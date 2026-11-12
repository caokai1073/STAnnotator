; Inno Setup 6 script of the STAnnotator installer (Windows 10 / 11, x64). packaging\windows\build.ps1 compiles it
; after PyInstaller:  ISCC /DAppVersion=<version> packaging\windows\installer.iss  ->  dist\STAnnotator-Setup-<version>.exe
;
; The user chooses the install folder (any drive; by default for this user only, no administrator rights needed) and,
; on a first installation, the data folder, written to %APPDATA%\STAnnotator\settings.json, which the app reads
; (config.settings_file). Uninstalling removes the program only: samples, annotations and the settings stay.

#ifndef AppVersion
  #define AppVersion "0.0.0"
#endif
#define AppName "STAnnotator"

[Setup]
AppId={{6F0C2B7E-5C1D-4A8E-9B3F-2D7A41C9E815}
AppName={#AppName}
AppVersion={#AppVersion}
AppVerName={#AppName} {#AppVersion}
AppPublisher=STAnnotator
AppPublisherURL=https://github.com/caokai1073/STAnnotator
DefaultDirName={autopf}\{#AppName}
DefaultGroupName={#AppName}
DisableDirPage=no
DisableProgramGroupPage=yes
PrivilegesRequired=lowest
PrivilegesRequiredOverridesAllowed=dialog
ArchitecturesAllowed=x64compatible
ArchitecturesInstallIn64BitMode=x64compatible
OutputDir=..\..\dist
OutputBaseFilename=STAnnotator-Setup-{#AppVersion}
SetupIconFile=STAnnotator.ico
UninstallDisplayIcon={app}\STAnnotator.exe
UninstallDisplayName={#AppName} {#AppVersion}
Compression=lzma2
SolidCompression=yes
WizardStyle=modern
CloseApplications=yes

[Tasks]
Name: "desktopicon"; Description: "Create a desktop shortcut"; GroupDescription: "Shortcuts:"

; An update first removes the program's libraries of the version it replaces: files the new version no longer has
; would stay otherwise (a library's renamed or moved module could then be loaded, and two versions' metadata give the
; bug report a wrong version). Not when the data folder was put inside it (DataInsideProgram).
[InstallDelete]
Type: filesandordirs; Name: "{app}\_internal"; Check: not DataInsideProgram

[Files]
Source: "..\..\dist\STAnnotator\*"; DestDir: "{app}"; Flags: ignoreversion recursesubdirs createallsubdirs

[Icons]
Name: "{group}\{#AppName}"; Filename: "{app}\STAnnotator.exe"
Name: "{group}\Uninstall {#AppName}"; Filename: "{uninstallexe}"
Name: "{autodesktop}\{#AppName}"; Filename: "{app}\STAnnotator.exe"; Tasks: desktopicon

[Run]
Filename: "{app}\STAnnotator.exe"; Description: "Start {#AppName}"; Flags: nowait postinstall skipifsilent

[Code]
var
  DataPage: TInputDirWizardPage;

function SettingsFile: String;
begin
  Result := ExpandConstant('{userappdata}\STAnnotator\settings.json');
end;

// The data folder of settings.json (JSON: backslashes doubled) inside the program's _internal folder: keep it
function DataInsideProgram: Boolean;
var
  Text: AnsiString;
  Inner: String;
begin
  Result := False;
  if LoadStringFromFile(SettingsFile, Text) then
  begin
    Inner := ExpandConstant('{app}') + '\_internal';
    StringChangeEx(Inner, '\', '\\', True);
    Result := Pos(Lowercase(Inner), Lowercase(String(Text))) > 0;
  end;
end;

procedure InitializeWizard;
begin
  DataPage := CreateInputDirPage(wpSelectDir, 'Data folder',
    'Where should STAnnotator keep its data?',
    'Imported samples, annotations and the images made for them go into this folder; a sample can take ' +
    'several GB, so choose a drive with room. It can be changed later in the app (? > Data folder > Change).',
    False, '');
  DataPage.Add('');
  DataPage.Values[0] := ExpandConstant('{%USERPROFILE}\stagent');
end;

function ShouldSkipPage(PageID: Integer): Boolean;
begin
  { an update keeps the data folder already chosen }
  Result := (PageID = DataPage.ID) and FileExists(SettingsFile);
end;

procedure CurStepChanged(CurStep: TSetupStep);
var
  Path: String;
  Lines: TArrayOfString;
begin
  if (CurStep = ssPostInstall) and not FileExists(SettingsFile) then
  begin
    Path := DataPage.Values[0];
    StringChangeEx(Path, '\', '\\', True);
    StringChangeEx(Path, '"', '\"', True);
    SetArrayLength(Lines, 3);
    Lines[0] := '{';
    Lines[1] := ' "data_folder": "' + Path + '"';
    Lines[2] := '}';
    ForceDirectories(ExtractFileDir(SettingsFile));
    ForceDirectories(DataPage.Values[0]);
    SaveStringsToUTF8File(SettingsFile, Lines, False);
  end;
end;
