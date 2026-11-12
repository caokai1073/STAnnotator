// STAnnotator bug reports: a Google Apps Script web app that saves each report the app uploads (Report bug... >
// Upload to the developers) into one Google Drive folder of the developer.
//
// Set up (once): create a Drive folder, put its id below (the part of the folder's URL after /folders/), paste this
// file into a new project at script.google.com, then Deploy > New deployment > Web app, "Execute as: Me",
// "Who has access: Anyone". Put the web app URL (https://script.google.com/macros/s/.../exec) in
// stagent/data/report_upload_url.txt before building the app (that file is not in git), so the built app uploads
// without asking; a tester running from source pastes it once. The URL is the only key: keep it out of public places.
// The receiver can only add files to the folder; a link that leaked is replaced by a new deployment.
// Automatic usage reports (the user allowed usage sharing) go into autosave/<installation id>/, made when missing;
// reports the user sends with Report a bug go into the folder itself. The screenshots a user adds to a report come
// as files of their own (one request each, "attachment": true), named <report name>-screenshot<i>.<ext>, next to it.
// The daily count ("checkin": true; every installation once a day: its random id, version and system) adds one row to
// the sheet STAnnotator-checkins in the folder, made by setup() (run it once in the editor: it asks for access to
// Sheets, which the count needs; then deploy a new version).

const FOLDER_ID = "PASTE_THE_FOLDER_ID_HERE";
const MAX_BYTES = 35 * 1024 * 1024;  // one request: Apps Script takes up to about 50 MB (the data is sent as base64);
                                     // with its screenshots apart, a report can be well over 100 MB in all
const ZONE = "America/New_York";     // the developers' day and clock (Boston): the count's days, its times, the caps
const MAX_PER_DAY = 100;             // reports a day, whoever sends them; screenshots are not counted
const MAX_CHECKINS_PER_DAY = 5000;  // rows a day in the count sheet, whoever sends them
const IMAGE_TYPES = ["image/png", "image/jpeg", "image/gif", "image/webp", "image/bmp", "image/tiff", "image/heic"];

function doPost(e) {
  const lock = LockService.getScriptLock();
  try {
    lock.waitLock(20000);
    const props = PropertiesService.getScriptProperties();
    const key = "count_" + Utilities.formatDate(new Date(), ZONE, "yyyy-MM-dd");
    const n = Number(props.getProperty(key) || 0);
    const body = JSON.parse(e.postData.contents);
    if (body.checkin) return checkin(body, props);
    const shot = !!body.attachment && !body.auto;
    if (n >= MAX_PER_DAY + (shot ? 1 : 0)) return reply({ error: "too many reports today: save the report to a file and send it" });
    const bytes = Utilities.base64Decode(body.data || "");
    if (!bytes.length) return reply({ error: shot ? "the screenshot is empty" : "the report is empty" });
    if (bytes.length > MAX_BYTES) return reply({ error: (shot ? "the screenshot" : "the report") + " is larger than 35 MB" });
    const name = String(body.name || "STAnnotator-report.zip").replace(/[^A-Za-z0-9_.-]/g, "_").slice(0, 120);
    const type = shot && IMAGE_TYPES.indexOf(body.type) >= 0 ? body.type : shot ? "application/octet-stream" : "application/zip";
    const file = folderFor(body).createFile(Utilities.newBlob(bytes, type, name));
    file.setDescription([`STAnnotator ${body.version || "?"} · ${body.system || ""}` + (body.install ? ` · installation ${body.install}` : ""),
                         body.auto ? "automatic usage report" : "",
                         shot ? `screenshot of ${String(body.report || "a report").slice(0, 120)}` : "", body.from || "",
                         String(body.note || "").slice(0, 2000)].filter(String).join("\n"));
    if (!shot) {
      Object.keys(props.getProperties()).filter((k) => k.startsWith("count_") && k !== key).forEach((k) => props.deleteProperty(k));
      props.setProperty(key, String(n + 1));
    }
    return reply({ ok: true, name: file.getName() });
  } catch (err) {
    return reply({ error: String(err) });
  } finally {
    lock.releaseLock();
  }
}

// The daily count: one row per installation and day (date, day the app counted, installation, version, system)
function checkin(body, props) {
  const key = "checkins_" + Utilities.formatDate(new Date(), ZONE, "yyyy-MM-dd");
  const n = Number(props.getProperty(key) || 0);
  if (n >= MAX_CHECKINS_PER_DAY) return reply({ error: "too many today" });
  const clean = (v, k) => String(v || "").replace(/[^A-Za-z0-9_.+ -]/g, "").slice(0, k);
  // the time as text in Boston time (a date cell shows in the sheet's own time zone); day: the app's day (Boston from
  // 0.2.57, UTC before); a row per day and installation, and another when it was updated that day
  sheet(props).appendRow([Utilities.formatDate(new Date(), ZONE, "yyyy-MM-dd HH:mm:ss"), clean(body.day, 10), clean(body.install, 40), clean(body.version, 20), clean(body.system, 20)]);
  Object.keys(props.getProperties()).filter((k) => k.startsWith("checkins_") && k !== key).forEach((k) => props.deleteProperty(k));
  props.setProperty(key, String(n + 1));
  return reply({ ok: true });
}

function sheet(props) {
  const id = props.getProperty("checkin_sheet");
  if (id) {
    try { return SpreadsheetApp.openById(id).getSheets()[0]; } catch (err) { /* deleted: made again */ }
  }
  const ss = SpreadsheetApp.create("STAnnotator-checkins");
  const sh = ss.getSheets()[0];
  sh.appendRow(["received (Boston)", "day", "installation", "version", "system"]);
  sh.setFrozenRows(1);
  DriveApp.getFileById(ss.getId()).moveTo(DriveApp.getFolderById(FOLDER_ID));
  props.setProperty("checkin_sheet", ss.getId());
  return sh;
}

// Run once in the editor after pasting this file: it asks for the access the receiver needs (Drive, Sheets) and makes
// the count sheet in the folder.
function setup() {
  const sh = sheet(PropertiesService.getScriptProperties());
  Logger.log("count sheet: " + sh.getParent().getUrl());
}

// Run once in the editor: the "day" of every row from the time it arrived, in Boston (ZONE). Apps before 0.2.57 sent
// their UTC day, so a count made in a Boston evening shows the next day. Also names the first column for Boston time.
// Rows whose arrival time is not a date (or a "yyyy-MM-dd ..." text) keep their day. Logs how many days it changed.
function bostonDays() {
  const sh = sheet(PropertiesService.getScriptProperties());
  const n = sh.getLastRow() - 1;
  sh.getRange(1, 1).setValue("received (Boston)");
  if (n < 1) return;
  const rows = sh.getRange(2, 1, n, 2).getValues();
  let changed = 0;
  const days = rows.map(([rec, day]) => {
    let d = String(day);
    if (rec instanceof Date && !isNaN(rec)) d = Utilities.formatDate(rec, ZONE, "yyyy-MM-dd");
    else if (/^\d{4}-\d{2}-\d{2}/.test(String(rec))) d = String(rec).slice(0, 10);
    else if (day instanceof Date) d = Utilities.formatDate(day, ZONE, "yyyy-MM-dd");
    if (d !== (day instanceof Date ? Utilities.formatDate(day, ZONE, "yyyy-MM-dd") : String(day))) changed++;
    return [d];
  });
  // (no number format: a table's typed column refuses one, and Apps Script raises it only when the changes are
  // written, past any try; a date-typed column shows the same day)
  sh.getRange(2, 2, n, 1).setValues(days);
  Logger.log(changed + " of " + n + " rows got another day");
}

function folderFor(body) {
  const root = DriveApp.getFolderById(FOLDER_ID);
  if (!body.auto) return root;
  const id = String(body.install || "").replace(/[^A-Za-z0-9_-]/g, "").slice(0, 40) || "unknown";
  return child(child(root, "autosave"), id);
}

function child(parent, name) {
  const found = parent.getFoldersByName(name);
  return found.hasNext() ? found.next() : parent.createFolder(name);
}

// The app's "Test" of the link
function doGet() {
  return reply({ ok: true, service: "STAnnotator reports" });
}

function reply(o) {
  return ContentService.createTextOutput(JSON.stringify(o)).setMimeType(ContentService.MimeType.JSON);
}
