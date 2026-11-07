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
// reports the user sends with Report a bug go into the folder itself.

const FOLDER_ID = "PASTE_THE_FOLDER_ID_HERE";
const MAX_BYTES = 35 * 1024 * 1024;  // Apps Script takes requests up to about 50 MB (the report is sent as base64)
const MAX_PER_DAY = 100;             // reports a day (UTC), whoever sends them

function doPost(e) {
  const lock = LockService.getScriptLock();
  try {
    lock.waitLock(20000);
    const props = PropertiesService.getScriptProperties();
    const key = "count_" + Utilities.formatDate(new Date(), "UTC", "yyyy-MM-dd");
    const n = Number(props.getProperty(key) || 0);
    if (n >= MAX_PER_DAY) return reply({ error: "too many reports today: save the report to a file and send it" });
    const body = JSON.parse(e.postData.contents);
    const bytes = Utilities.base64Decode(body.data || "");
    if (!bytes.length) return reply({ error: "the report is empty" });
    if (bytes.length > MAX_BYTES) return reply({ error: "the report is larger than 35 MB" });
    const name = String(body.name || "STAnnotator-report.zip").replace(/[^A-Za-z0-9_.-]/g, "_").slice(0, 120);
    const file = folderFor(body).createFile(Utilities.newBlob(bytes, "application/zip", name));
    file.setDescription([`STAnnotator ${body.version || "?"} · ${body.system || ""}` + (body.install ? ` · installation ${body.install}` : ""),
                         body.auto ? "automatic usage report" : "", body.from || "",
                         String(body.note || "").slice(0, 2000)].filter(String).join("\n"));
    Object.keys(props.getProperties()).filter((k) => k.startsWith("count_") && k !== key).forEach((k) => props.deleteProperty(k));
    props.setProperty(key, String(n + 1));
    return reply({ ok: true, name: file.getName() });
  } catch (err) {
    return reply({ error: String(err) });
  } finally {
    lock.releaseLock();
  }
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
