// STAnnotator bug reports: a Google Apps Script web app that saves each report the app uploads (Report bug... >
// Upload to the developers) into one Google Drive folder of the developer.
//
// Set up (once): create a Drive folder, put its id below (the part of the folder's URL after /folders/), paste this
// file into a new project at script.google.com, then Deploy > New deployment > Web app, "Execute as: Me",
// "Who has access: Anyone". Give the web app URL (https://script.google.com/macros/s/.../exec) to the testers; they
// paste it once into the app. The URL is the only key: keep it out of public places.

const FOLDER_ID = "PASTE_THE_FOLDER_ID_HERE";
const MAX_BYTES = 40 * 1024 * 1024;  // Apps Script takes requests up to about 50 MB (the report is sent as base64)

function doPost(e) {
  try {
    const body = JSON.parse(e.postData.contents);
    const bytes = Utilities.base64Decode(body.data || "");
    if (!bytes.length) return reply({ error: "the report is empty" });
    if (bytes.length > MAX_BYTES) return reply({ error: "the report is larger than 40 MB" });
    const name = String(body.name || "STAnnotator-report.zip").replace(/[^A-Za-z0-9_.-]/g, "_").slice(0, 120);
    const file = DriveApp.getFolderById(FOLDER_ID).createFile(Utilities.newBlob(bytes, "application/zip", name));
    file.setDescription([`STAnnotator ${body.version || "?"} · ${body.system || ""}`, body.from || "",
                         String(body.note || "").slice(0, 2000)].filter(String).join("\n"));
    return reply({ ok: true, name: file.getName() });
  } catch (err) {
    return reply({ error: String(err) });
  }
}

// The app's "Test" of the link
function doGet() {
  return reply({ ok: true, service: "STAnnotator reports" });
}

function reply(o) {
  return ContentService.createTextOutput(JSON.stringify(o)).setMimeType(ContentService.MimeType.JSON);
}
