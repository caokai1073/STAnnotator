/* STAnnotator interactive tutorial: the player.
 *
 * Stage: the replayed Claude Code terminal on the left, the real app (stagent/static/app.html, unchanged) on the right
 * in a same-origin iframe. The app's requests are answered by shim.js from a recorded session (data/manifest.json),
 * per scene: a scene is a named checkpoint of that session (after the import, after the alignment, ...). A step can
 * switch the scene, so the app shows what the real app showed at that moment.
 *
 * Script: chapters.js (window.TUTORIAL). Actions run on the app's own page: a virtual cursor moves to the element and
 * clicks it for real, so the app responds as it does for a user. Pausing removes the shield: the app is the viewer's.
 */
(() => {
  "use strict";
  const T = window.TUTORIAL;
  const $ = (s, r = document) => r.querySelector(s);
  const el = (tag, cls, html) => { const e = document.createElement(tag); if (cls) e.className = cls; if (html != null) e.innerHTML = html; return e; };
  const esc = (s) => String(s).replace(/[&<>]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;" }[c]));
  const rich = (s) => esc(s).replace(/\*\*(.+?)\*\*/g, "<b>$1</b>").replace(/`([^`]+)`/g, "<code>$1</code>");
  // the agent's markdown: tables become tables (as Claude Code draws them), the rest stays text
  const richMd = (s) => {
    const out = [], lines = String(s).split("\n");
    for (let i = 0; i < lines.length; i++) {
      if (/^\s*\|.*\|\s*$/.test(lines[i])) {
        const rows = [];
        while (i < lines.length && /^\s*\|.*\|\s*$/.test(lines[i])) { rows.push(lines[i]); i++; }
        i--;
        const cells = rows.filter((r) => !/^\s*\|[\s:|-]+\|\s*$/.test(r)).map((r) => r.trim().slice(1, -1).split("|").map((c) => rich(c.trim())));
        out.push('<table class="t-table">' + cells.map((r, k) => "<tr>" + r.map((c) => (k ? "<td>" : "<th>") + c + (k ? "</td>" : "</th>")).join("") + "</tr>").join("") + "</table>");
      } else out.push(rich(lines[i]) + "\n");
    }
    return out.join("");
  };
  const fmt = (ms) => { const s = Math.max(0, Math.round(ms / 1000)); return `${Math.floor(s / 60)}:${String(s % 60).padStart(2, "0")}`; };
  const CANCEL = Symbol("cancel");
  const SPEEDS = [0.75, 1, 1.25, 1.5, 2];
  const DATA = window.TUTORIAL_DATA || "data/";  // the recorded session: next to the page, or the media repository
  const BLANK = "data:image/gif;base64,R0lGODlhAQABAAAAACH5BAEKAAEALAAAAAABAAEAAAICTAEAOw==";

  /* ---------------- state ---------------- */
  const S = { ci: 0, si: 0, playing: false, run: 0, speed: 1, scene: 0, data: null, agent: {}, dirty: false,
              stepStart: 0, stepElapsed: 0, ready: false };
  const frame = $("#frame");
  const W = () => frame.contentWindow;
  const D = () => frame.contentDocument;

  /* ---------------- replay of the recorded answers (used by shim.js inside the app page) ---------------- */
  const counters = new Map();
  const keyOf = (method, url) => {
    const u = new URL(url, "http://x");
    const q = [...u.searchParams].filter(([k]) => k !== "t" && k !== "_")
      .sort((a, b) => (a[0] < b[0] ? -1 : a[0] > b[0] ? 1 : a[1] < b[1] ? -1 : a[1] > b[1] ? 1 : 0));
    let path = u.pathname;
    try { path = decodeURIComponent(path); } catch (e) { /* keep it */ }  // the server records decoded paths
    return method + " " + path + (q.length ? "?" + q.map(([k, v]) => k + "=" + v).join("&") : "");
  };
  function match(method, url, body) {
    const key = keyOf(method, url), E = S.data.entries, names = S.data.scenes;
    for (let i = S.scene; i >= 0; i--) {
      const list = E[names[i]] && E[names[i]][key];
      if (!list) continue;
      let ans;
      if (method !== "GET") {
        ans = list.find((a) => a[3] === body) || list[Math.min(counters.get(i + key) || 0, list.length - 1)];
      } else if (/^GET \/api\/job\//.test(key)) {  // a job's progress: in the recorded order
        const c = counters.get(i + key) || 0;
        ans = list[Math.min(c, list.length - 1)];
      } else ans = list[list.length - 1];  // anything else: what it was at the end of the scene
      counters.set(i + key, (counters.get(i + key) || 0) + 1);
      return { status: ans[0], ctype: ans[1], href: new URL(DATA + "b/" + ans[2], location.href).href };
    }
    const bare = key.split("?")[0];  // H&E tiles are stored without their revision: the same image file
    const rev = /\/he\/\d+\//.test(bare) && new URL(url, "http://x").searchParams.get("r");
    const blank = rev && (S.data.blank_tiles || {})[rev];  // an image waiting for its alignment: a blank canvas, not the H&E
    if (blank) return { status: 200, ctype: "image/jpeg", href: new URL(DATA + "b/" + blank, location.href).href };
    for (const g of S.data.global || []) {  // answers for every scene: tiles, gene colors, the stand-ins for Terminal
      const list = E[g] && (E[g][key] || (/\/he\/\d+\//.test(bare) && E[g][bare]));
      if (list) { const a = list[list.length - 1]; return { status: a[0], ctype: a[1], href: new URL(DATA + "b/" + a[2], location.href).href }; }
    }
    for (let i = S.scene + 1; i < names.length; i++) {  // recorded only later (option lists and the like): the nearest one
      const list = E[names[i]] && E[names[i]][key];
      if (list) { const a = list[0]; return { status: a[0], ctype: a[1], href: new URL(DATA + "b/" + a[2], location.href).href }; }
    }
    if (S.debug) console.warn("[tutorial] not recorded:", key);
    if (S.check) (S.misses = S.misses || []).push({ ch: CH[S.ci].id, step: S.si, key });
    return null;
  }
  function setScene(name) {
    const i = S.data.scenes.indexOf(name);
    if (i < 0) { console.warn("[tutorial] unknown scene", name); return; }
    if (i !== S.scene) { S.scene = i; counters.clear(); }
  }
  window.__TUT = {
    match, blank: BLANK,
    toast: (m) => toast(m),
    download: () => toast("Downloads need the installed app: this demo replays a recorded session."),
    attach: (w) => { w.addEventListener("pointerdown", markDirty, true); w.addEventListener("keydown", markDirty, true); },
    log: (m, u, hit) => { if (S.debug) console.log("[tutorial]", hit ? "hit " : "MISS", m, u); },
  };
  function markDirty(e) { if (!S.playing && e.isTrusted) S.dirty = true; }

  /* ---------------- timing ---------------- */
  const readMs = (s) => Math.min(15000, Math.max(3200, 900 + 52 * String(s || "").replace(/\*\*/g, "").length));
  function check(run) { if (run !== S.run) throw CANCEL; }
  async function sleep(ms, run, instant) {  // counts down only while playing; speed applies
    if (instant) return check(run);
    let left = ms, last = performance.now();
    while (left > 0) {
      await new Promise((r) => setTimeout(r, 40));
      check(run);
      const now = performance.now();
      if (S.playing) left -= (now - last) * S.speed;
      last = now;
    }
  }
  async function until(fn, run, timeout = 8000) {
    const t0 = performance.now();
    while (performance.now() - t0 < timeout) {
      check(run);
      try { if (fn()) return true; } catch (e) { /* not there yet */ }
      await new Promise((r) => setTimeout(r, 60));
    }
    return false;
  }

  /* ---------------- targets in the app page ---------------- */
  function find(spec) {
    const d = D();
    if (!d) return null;
    if (typeof spec === "function") return spec(d, W());
    if (typeof spec === "string") return visible(d.querySelectorAll(spec))[0] || null;
    let list = spec.id ? [d.getElementById(spec.id)].filter(Boolean) : [...d.querySelectorAll(spec.css || "button,[role=button],a,label,select,input,textarea,summary,div,span")];
    if (spec.title) list = list.filter((e) => (e.getAttribute("title") || "").startsWith(spec.title));
    if (spec.text != null) {
      const t = spec.text;
      list = list.filter((e) => { const s = (e.textContent || "").trim().replace(/\s+/g, " "); return t instanceof RegExp ? t.test(s) : s === t; });
      if (!spec.css) list = list.filter((e) => ![...e.children].some((c) => list.includes(c)));  // innermost
    }
    if (spec.ph) list = list.filter((e) => (e.getAttribute("placeholder") || "").startsWith(spec.ph));
    list = visible(list);
    let e = list[spec.nth || 0] || null;
    if (e && spec.up) e = e.closest(spec.up);
    return e;
  }
  function visible(list) { return [...list].filter((e) => e && e.getClientRects().length && getComputedStyle(e).visibility !== "hidden"); }
  function rectOf(target) {  // in stage coordinates of the app window
    if (target && target.box) return target.box;
    const e = target instanceof Element ? target : find(target);
    if (!e) return null;
    const r = e.getBoundingClientRect();
    return { x: r.left, y: r.top, w: r.width, h: r.height };
  }

  /* ---------------- effects: cursor, spotlight, ripple ---------------- */
  const cur = $("#cursor"), spot = $("#spot"), spotLabel = $("#spotLabel"), ripple = $("#ripple");
  let curPos = { x: 700, y: 520 };
  async function moveTo(x, y, run, instant) {
    cur.classList.add("on");
    const dist = Math.hypot(x - curPos.x, y - curPos.y);
    const ms = instant ? 0 : Math.min(1100, 260 + dist * 0.9) / S.speed;
    cur.style.transition = instant ? "none" : `left ${ms}ms cubic-bezier(.4,.1,.2,1), top ${ms}ms cubic-bezier(.4,.1,.2,1)`;
    cur.style.left = x + "px"; cur.style.top = y + "px";
    curPos = { x, y };
    if (!instant) await sleep(ms * S.speed + 80, run);
  }
  function showSpot(target, label) {
    const r = target && rectOf(target);
    if (!r) { hideSpot(); return; }
    const pad = 5;
    Object.assign(spot.style, { left: r.x - pad + "px", top: r.y - pad + "px", width: r.w + 2 * pad + "px", height: r.h + 2 * pad + "px" });
    spot.classList.add("on");
    if (label) {
      spotLabel.textContent = label;
      spotLabel.classList.add("on");
      const appW = frame.clientWidth, lw = spotLabel.offsetWidth || 160;
      let lx = r.x + r.w + 14, ly = r.y + r.h / 2 - 13;
      if (lx + lw > appW - 8) { lx = Math.max(8, r.x + r.w / 2 - lw / 2); ly = r.y + r.h + 12; if (ly > frame.clientHeight - 40) ly = r.y - 38; }
      spotLabel.style.left = lx + "px"; spotLabel.style.top = ly + "px";
    } else spotLabel.classList.remove("on");
  }
  function hideSpot() { spot.classList.remove("on"); spotLabel.classList.remove("on"); }
  function clickFx() { ripple.style.left = curPos.x + "px"; ripple.style.top = curPos.y + "px"; ripple.classList.remove("go"); void ripple.offsetWidth; ripple.classList.add("go"); }

  /* ---------------- actions ---------------- */
  function realClick(e) {
    const r = e.getBoundingClientRect(), o = { bubbles: true, cancelable: true, view: W(), clientX: r.left + r.width / 2, clientY: r.top + r.height / 2 };
    e.dispatchEvent(new (W().PointerEvent || W().MouseEvent)("pointerdown", o));
    e.dispatchEvent(new (W().MouseEvent)("mousedown", o));
    e.dispatchEvent(new (W().PointerEvent || W().MouseEvent)("pointerup", o));
    e.dispatchEvent(new (W().MouseEvent)("mouseup", o));
    e.click();
  }
  async function pointAt(target, run, instant) {
    const e = target instanceof Element ? target : find(target);
    if (!e) throw new Error("target not found: " + JSON.stringify(target, (k, v) => (v instanceof RegExp ? String(v) : v)));
    if (!instant) e.scrollIntoView && e.scrollIntoView({ block: "nearest" });
    const r = e.getBoundingClientRect();
    await moveTo(r.left + Math.min(r.width / 2, 40 + r.width * 0.25), r.top + r.height / 2, run, instant);
    return e;
  }
  const ACT = {
    async click(run, instant, target) {
      const e = await pointAt(target, run, instant);
      if (!instant) { clickFx(); await sleep(160, run); }
      realClick(e);
      if (!instant) await sleep(420, run);
    },
    async select(run, instant, target, value) {
      const e = await pointAt(target, run, instant);
      if (!instant) { clickFx(); await fakeMenu(e, value, run); }
      const opt = [...e.options].find((o) => o.value === value || o.textContent.trim() === value);
      e.value = opt ? opt.value : value;
      e.dispatchEvent(new (W().Event)("input", { bubbles: true }));
      e.dispatchEvent(new (W().Event)("change", { bubbles: true }));
      if (!instant) await sleep(500, run);
    },
    async type(run, instant, target, text) {
      const e = await pointAt(target, run, instant);
      if (!instant) { clickFx(); await sleep(150, run); }
      e.focus && e.focus();
      if (instant) { e.value = text; }
      else {
        e.value = "";
        const per = Math.max(18, Math.min(55, 2400 / text.length));
        for (const ch of text) { e.value += ch; e.dispatchEvent(new (W().Event)("input", { bubbles: true })); await sleep(per, run); }
      }
      e.dispatchEvent(new (W().Event)("input", { bubbles: true }));
      e.dispatchEvent(new (W().Event)("change", { bubbles: true }));
      if (!instant) await sleep(300, run);
    },
    async hover(run, instant, target) { await pointAt(target, run, instant); },
    async move(run, instant, x, y) { await moveTo(x, y, run, instant); },
    async wait(run, instant, ms) { await sleep(ms, run, instant); },
    async scene(run, instant, name) { setScene(name); },
    async call(run, instant, fn, ...args) { const f = W()[fn]; if (typeof f === "function") await Promise.resolve(f(...args)).catch(() => {}); },
    async eval(run, instant, code) {
      try { await Promise.resolve(W().eval(code)); } catch (e) { console.warn("[tutorial] eval", e); }
    },
    async until(run, instant, code, timeout) { await until(() => W().eval(code), run, timeout || 10000); },
    async spot(run, instant, target, label) { if (!instant) showSpot(target, label); },
    async unspot() { hideSpot(); },
    async zoom(run, instant, box, ms) {  // box in image pixels of the H&E: {x, y, w, h}, or "home"
      const v = W().eval("typeof viewer !== 'undefined' ? viewer : null");
      if (!v) return;
      if (box === "home") v.viewport.goHome(instant);
      else {
        const it = v.world.getItemAt(0);
        if (!it) return;
        const r = it.imageToViewportRectangle(box.x, box.y, box.w, box.h);
        v.viewport.fitBounds(r, instant);
      }
      if (!instant) await sleep(ms || 1600, run);
    },
    async agent(run, instant, id, opt) { await playAgent(id, opt || {}, run, instant); },
    async editDrag(run, instant, moves) {  // move points of the outline being edited: [[x, y, dy] in H&E pixels]
      const w = W();
      for (const [x, y, dy] of moves) {
        check(run);
        const e = w.eval("st.editing");
        if (!e) return;
        let best = 0, bd = Infinity;
        e.pts.forEach((q, i) => { const d = Math.hypot(q[0] - x, q[1] - y); if (d < bd) { bd = d; best = i; } });
        if (instant) { e.pts[best] = [e.pts[best][0], e.pts[best][1] + dy]; e.dirty = true; w.eval("drawRegions()"); continue; }
        const a = w.toScreen(e.pts[best][0], e.pts[best][1]), b = w.toScreen(e.pts[best][0], e.pts[best][1] + dy);
        const osd = D().getElementById("osd").getBoundingClientRect();
        const ax = osd.left + a.x, ay = osd.top + a.y, bx = osd.left + b.x, by = osd.top + b.y;
        await moveTo(ax, ay, run);
        const hdl = [...D().querySelectorAll("#regsvg circle.hdl:not(.mid)")].sort((p, q) => {
          const rp = p.getBoundingClientRect(), rq = q.getBoundingClientRect();
          return Math.hypot(rp.left + rp.width / 2 - ax, rp.top + rp.height / 2 - ay) - Math.hypot(rq.left + rq.width / 2 - ax, rq.top + rq.height / 2 - ay);
        })[0];
        const ev = (t, px, py, target) => target.dispatchEvent(new (w.PointerEvent)(t, { bubbles: true, cancelable: true, clientX: px, clientY: py, pointerId: 1, pointerType: "mouse", isPrimary: true, buttons: 1 }));
        if (!hdl) continue;
        clickFx(); ev("pointerdown", ax, ay, hdl);
        cur.style.transition = "none";
        for (let k = 1; k <= 14; k++) {
          const px = ax + (bx - ax) * k / 14, py = ay + (by - ay) * k / 14;
          cur.style.left = px + "px"; cur.style.top = py + "px"; curPos = { x: px, y: py };
          ev("pointermove", px, py, w); await sleep(45, run);
        }
        ev("pointerup", bx, by, w); await sleep(350, run);
      }
    },
    async picker(run, instant, opt) {  // the desktop app's file dialog (a native window): drawn here, for the viewer to see the file
      if (instant) return;
      const box = el("div", "picker");
      box.innerHTML = `<div class="ph"><span class="dots"><i></i><i></i><i></i></span><b>${esc(opt.title || "Open")}</b></div>
        <div class="pp"><span>⌂ ${esc(opt.folder)}</span></div><div class="pl"></div>
        <div class="pf"><span>${esc(opt.kinds || "")}</span><button class="pc">Cancel</button><button class="po">Open</button></div>`;
      const list = box.querySelector(".pl");
      let pick = null;
      for (const f of opt.files) {
        const r = el("div", "pr" + (f === opt.pick ? "" : /\.(geojson|json)$/.test(f) ? "" : " dim"), `<span class="ic">${/\.\w+$/.test(f) ? "📄" : "📁"}</span>${esc(f)}`);
        if (f === opt.pick) pick = r;
        list.append(r);
      }
      $("#fx").append(box);
      try {
        await sleep(700, run);
        if (pick) {
          const fr = $("#fx").getBoundingClientRect(), k = fr.width / $("#fx").offsetWidth, pr = pick.getBoundingClientRect();
          await moveTo((pr.left - fr.left) / k + 60, (pr.top - fr.top) / k + pr.height / k / 2, run);
          clickFx(); pick.classList.add("sel"); await sleep(opt.hold || 3500, run);
          const ob = box.querySelector(".po").getBoundingClientRect();
          await moveTo((ob.left - fr.left) / k + ob.width / k / 2, (ob.top - fr.top) / k + ob.height / k / 2, run);
          clickFx(); await sleep(300, run);
        }
      } finally { box.remove(); }
    },

    async term(run, instant, what) { if (what === "clear") termReset(); },
  };
  async function fakeMenu(sel, value, run) {  // a native <select> cannot be opened by script: draw its list
    const r = sel.getBoundingClientRect(), fx = $("#fx");
    const m = el("div");
    Object.assign(m.style, { position: "absolute", left: r.left + "px", top: r.bottom + 2 + "px", minWidth: r.width + "px", background: "#fff",
      border: "1px solid #d0d4db", borderRadius: "8px", boxShadow: "0 8px 24px rgba(0,0,0,.18)", padding: "4px", font: "13px -apple-system, sans-serif", zIndex: 5, maxHeight: "300px", overflow: "hidden" });
    let pick = null;
    for (const o of [...sel.options].slice(0, 14)) {
      const it = el("div", "", esc(o.textContent));
      Object.assign(it.style, { padding: "4px 10px", borderRadius: "5px", whiteSpace: "nowrap", color: "#1d2433" });
      if (o.value === value || o.textContent.trim() === value) pick = it;
      m.append(it);
    }
    fx.append(m);
    try {
      await sleep(450, run);
      if (pick) {
        const pr = pick.getBoundingClientRect(), fr = fx.getBoundingClientRect(), k = fr.width / fx.offsetWidth;
        await moveTo((pr.left - fr.left) / k + 30, (pr.top - fr.top) / k + pr.height / k / 2, run);
        pick.style.background = "#fff2dc"; await sleep(260, run); clickFx(); await sleep(200, run);
      }
    } finally { m.remove(); }
  }

  /* ---------------- the replayed terminal ---------------- */
  const tbody = $("#tbody");
  let tscroll, tprompt;
  function termReset() {
    tbody.textContent = "";
    tscroll = el("div", "scroll");
    tscroll.append(el("div", "t-welcome", "<b>✻ Claude Code</b><div>cwd: ~/stagent/claude_workspace · stagent tools connected</div>"));
    tprompt = el("div", "t-prompt", '<span class="typed"></span><span class="caret"></span>');
    tbody.append(tscroll, tprompt);
  }
  const shortArgs = (o) => Object.entries(o || {}).map(([k, v]) => `${k}: ${typeof v === "string" ? JSON.stringify(v.length > 60 ? v.slice(0, 57) + "…" : v) : JSON.stringify(v).slice(0, 60)}`).join(", ");
  function outLines(s) {  // as Claude Code shows a tool's answer: its first lines, then how many more
    s = String(s || "").trim();
    if (!s) return "(no output)";
    let j = null;
    try { j = JSON.parse(s); if (j && typeof j.result === "string") { try { j = JSON.parse(j.result); } catch (e) { s = j.result; j = null; } } } catch (e) { /* plain text */ }
    const lines = (j ? JSON.stringify(j, null, 2) : s).split("\n").map((l) => (l.length > 62 ? l.slice(0, 60) + "…" : l));
    const n = s.includes("[Image]") ? 2 : 3;
    return lines.slice(0, n).join("\n") + (lines.length > n ? `\n… +${lines.length - n} lines (ctrl+o to expand)` : "") + (s.includes("[Image]") ? "\n[Image]" : "");
  }
  function add(node) {
    tscroll.append(node);
    while (tscroll.children.length > 60) tscroll.children[1].remove();
    tscroll.scrollTop = tscroll.scrollHeight;
  }
  function agentMs(ev, opt) {  // estimated playing time of one conversation at 1x
    let ms = 0;
    for (const e of ev) {
      if (e.k === "user") ms += opt.fromApp ? 900 : Math.min(3500, 28 * e.text.length) + 400;
      else if (e.k === "tool") ms += 1500;
      else if (e.k === "text") ms += Math.min(2600, 10 * e.text.length) + 900;
    }
    const last = [...ev].reverse().find((e) => e.k === "text");
    return ms + (last ? Math.min(9000, 22 * last.text.length) : 0);
  }
  async function playAgent(id, opt, run, instant) {
    const ev = S.agent[id];
    if (!ev) { console.warn("[tutorial] no conversation", id); return; }
    const scenes = opt.scenes || {};
    let spin = null;
    const spinOn = (label) => { if (instant) return; if (!spin) { spin = el("div", "t-spin"); add(spin); } spin.textContent = "✻ " + (label || "Working…"); tscroll.append(spin); };
    const spinOff = () => { if (spin) { spin.remove(); spin = null; } };
    for (const e of ev) {
      check(run);
      if (e.k === "user") {
        const typed = tprompt.querySelector(".typed");
        if (!instant && !opt.fromApp) {
          const per = Math.max(12, Math.min(40, 3500 / e.text.length));
          for (const ch of e.text) { typed.textContent += ch; await sleep(per, run); }
          await sleep(350, run);
        }
        typed.textContent = "";
        add(el("div", "t-user" + (opt.fromApp ? " from-app" : ""), esc(opt.show || e.text)));
        spinOn("Thinking…");
        await sleep(opt.fromApp ? 900 : 600, run, instant);
      } else if (e.k === "tool") {
        const nm = e.name.replace(/^mcp__stagent__/, "");
        const line = el("div", "t-line t-tool" + (e.error ? " err" : ""), `<span class="nm">stagent · ${esc(nm)}</span><span class="args">(${esc(shortArgs(e.input))})</span>`);
        add(line); spinOn(); await sleep(750, run, instant);
        add(el("div", "t-out", esc(outLines(e.out)))); spinOn(); tscroll.scrollTop = tscroll.scrollHeight;
        if (scenes[nm]) setScene(scenes[nm]);
        await sleep(750, run, instant);
      } else if (e.k === "text") {
        spinOff();
        const node = el("div", "t-line t-text"); add(node);
        if (instant) node.innerHTML = richMd(e.text);
        else {
          const words = e.text.split(/(\s+)/), per = Math.max(8, Math.min(40, 2600 / words.length));
          let acc = "";
          for (const w of words) { acc += w; node.innerHTML = rich(acc); tscroll.scrollTop = tscroll.scrollHeight; await sleep(per, run); }
          node.innerHTML = richMd(e.text); tscroll.scrollTop = tscroll.scrollHeight;
        }
        if (e !== ev[ev.length - 2]) { spinOn(); await sleep(900, run, instant); }
      } else if (e.k === "done") {
        spinOff();
        add(el("div", "t-done", `✻ Worked for ${fmt(e.ms || 0)}`));
        if (scenes.$done) setScene(scenes.$done);
        const last = [...ev].reverse().find((x) => x.k === "text");
        await sleep(last ? Math.min(9000, 22 * last.text.length) : 1500, run, instant);
      }
    }
    spinOff();
  }

  /* ---------------- chapters, steps, timeline ---------------- */
  const CH = T.chapters;
  function stepMs(st) {
    let ms = readMs(st.say);
    let act = 0;
    for (const a of st.do || []) {
      if (a[0] === "click" || a[0] === "hover") act += 1300;
      else if (a[0] === "select") act += 2400;
      else if (a[0] === "type") act += 900 + Math.min(2400, 30 * String(a[2]).length);
      else if (a[0] === "wait") act += a[1];
      else if (a[0] === "zoom") act += a[2] || 1600;
      else if (a[0] === "agent") act += S.agent[a[1]] ? agentMs(S.agent[a[1]], a[2] || {}) : 0;
    }
    return st.ms || Math.max(ms, act + 1200);
  }
  let TL = [];  // [{ci, si, t0, ms}]
  function buildTimeline() {
    TL = []; let t = 0;
    CH.forEach((c, ci) => { c.t0 = t; c.steps.forEach((st, si) => { const ms = stepMs(st); TL.push({ ci, si, t0: t, ms }); t += ms; }); c.ms = t - c.t0; });
    S.total = t;
  }
  const tl = (ci, si) => TL.find((x) => x.ci === ci && x.si === si);

  /* ---------------- UI ---------------- */
  function buildNav() {
    const nav = $("#nav"); nav.textContent = "";
    let grp = null;
    CH.forEach((c, ci) => {
      if (c.group !== grp) { grp = c.group; nav.append(el("div", "grp", esc(grp))); }
      const b = el("button", "", `<span class="n">${ci + 1}</span><span>${esc(c.title)}</span><span class="d">${fmt(c.ms)}</span>`);
      b.onclick = () => jump(ci, 0, true);
      nav.append(b);
    });
    const segs = $("#segs"); segs.textContent = "";
    CH.forEach((c) => { const d = el("div"); d.style.flex = `${c.ms} 0 0`; segs.append(d); });
  }
  function renderNotes() {
    const c = CH[S.ci], st = c.steps[S.si];
    const n = $("#notes"); n.textContent = "";
    n.append(el("span", "chip", `Chapter ${S.ci + 1} of ${CH.length} · ${fmt(c.ms)}`));
    n.append(el("h2", "", esc(c.title)));
    n.append(el("p", "sum", rich(c.summary || "")));
    const card = el("div", "card", `<div class="k">Step ${S.si + 1} / ${c.steps.length}</div><h3>${esc(st.title)}</h3>`);
    if (st.pts && st.pts.length) card.append(el("ul", "", st.pts.map((p) => `<li>${rich(p)}</li>`).join("")));
    if (st.tip) card.append(el("div", "tip", "Tip · " + rich(st.tip)));
    n.append(card);
    const steps = el("div", "steps", '<div class="k">Steps in this chapter (click to jump)</div>');
    c.steps.forEach((s, si) => {
      const b = el("button", si < S.si ? "done" : si === S.si ? "cur" : "", `<span class="i"><span>${si + 1}</span></span><span>${esc(s.title)}</span>`);
      b.onclick = () => jump(S.ci, si, true);
      steps.append(b);
    });
    n.append(steps);
    [...$("#nav").querySelectorAll("button")].forEach((b, i) => { b.classList.toggle("cur", i === S.ci); b.classList.toggle("done", i < S.ci); });
    $("#caption").innerHTML = rich(st.say || "");
    $("#where").innerHTML = `${esc(c.title)} <span>· ${esc(st.title)}</span>`;
  }
  function progress() {
    const x = tl(S.ci, S.si); if (!x) return;
    const t = x.t0 + Math.min(x.ms, S.stepElapsed + (S.playing ? (performance.now() - S.stepStart) * S.speed : 0));
    const f = t / S.total;
    $("#fill").style.width = f * 100 + "%"; $("#knob").style.left = f * 100 + "%";
    $("#time").textContent = `${fmt(t)} / ${fmt(S.total)}`;
  }
  (function tick() { progress(); requestAnimationFrame(tick); })();
  // embedded in the website (?embed): the page's chapter buttons jump here, and show where the tutorial is
  const host = () => (S.embed && window.parent !== window ? window.parent : null);
  let lastPost = 0;
  function postState(force) {
    const h = host(), x = tl(S.ci, S.si);
    if (!h || !x || (!force && performance.now() - lastPost < 250)) return;
    lastPost = performance.now();
    const c = CH[S.ci], t = x.t0 + Math.min(x.ms, S.stepElapsed + (S.playing ? (performance.now() - S.stepStart) * S.speed : 0));
    try { h.postMessage({ source: "stannotator-tutorial", type: "state", chapter: c.id, progress: Math.max(0, Math.min(1, (t - c.t0) / c.ms)), playing: S.playing }, location.origin); } catch (e) { /* another site */ }
  }
  setInterval(() => postState(false), 300);
  window.addEventListener("message", (e) => {
    if (e.origin !== location.origin || !e.data || e.data.source !== "stannotator-site") return;
    if (e.data.type === "chapter") { const ci = CH.findIndex((c) => c.id === e.data.chapter); if (ci >= 0) jump(ci, 0, true); }
  });

  function setPlaying(p) {
    if (p && S.dirty) { S.dirty = false; S.playing = true; jump(S.ci, S.si, false); return; }  // the viewer changed the app: start this step again
    if (p) $("#app").classList.add("started");
    if (p === S.playing) return;
    if (p) S.stepStart = performance.now(); else S.stepElapsed += (performance.now() - S.stepStart) * S.speed;
    S.playing = p;
    $("#app").classList.toggle("paused", !p);
    $("#bPlay").textContent = p ? "❚❚" : "▶";
    if (!p) { hideSpot(); cur.classList.remove("on"); }
  }

  /* ---------------- loading a chapter, running steps ---------------- */
  async function loadApp(c, run) {
    $("#loading").classList.add("on");
    setScene(c.scene);
    counters.clear();
    S.ready = false;
    await new Promise((res) => { frame.onload = res; frame.src = "app/index.html?tut=" + Date.now(); });
    check(run);
    const G = (expr) => W().eval(expr);  // the app's top-level const (st, viewer) are not window properties
    await until(() => W().document.readyState === "complete" && G("typeof st === 'object' && Array.isArray(st.samples || [])"), run, 10000);
    if (c.sample) await until(() => G("st.sample && st.sample.sample_id"), run, 10000);
    await new Promise((r) => setTimeout(r, 400));
  }
  async function runActs(acts, run, instant) {
    for (const a of acts || []) {
      check(run);
      try { await ACT[a[0]](run, instant, ...a.slice(1)); }
      catch (e) {
        if (e === CANCEL) throw e;
        console.warn("[tutorial] step action failed", a, e);
        (S.fails = S.fails || []).push({ ch: CH[S.ci].id, step: S.si, act: a[0], err: String(e.message || e).slice(0, 160) });
      }
    }
  }
  async function jump(ci, si, autoplay) {
    const run = ++S.run;
    S.ci = ci; S.si = si; S.stepElapsed = 0; S.stepStart = performance.now();
    const keepPlaying = autoplay || S.playing;
    hideSpot(); cur.classList.remove("on");
    renderNotes();
    try {
      const c = CH[ci];
      termReset();
      // the agent's pane only in chapters where the agent works: elsewhere the app has the whole stage
      $("#stage").classList.toggle("noTerm", !c.steps.some((st) => (st.do || []).some((a) => a[0] === "agent")));
      fitApp();
      await loadApp(c, run);
      await runActs(c.pre, run, true);
      $("#loading").firstChild.textContent = si ? "Preparing this step…" : "Loading…";
      for (let k = 0; k < si; k++) { const st = c.steps[k]; if (st.scene) setScene(st.scene); await runActs(st.do, run, true); }
      $("#loading").classList.remove("on");
      if (keepPlaying) setPlaying(true);
      await play(run);
    } catch (e) { if (e !== CANCEL) console.error(e); }
  }
  async function play(run) {
    for (;;) {
      const c = CH[S.ci], st = c.steps[S.si];
      renderNotes();
      if (S.onStep) S.onStep();
      S.stepElapsed = 0; S.stepStart = performance.now();
      while (!S.playing) { await new Promise((r) => setTimeout(r, 60)); check(run); }
      if (st.scene) setScene(st.scene);
      const t0 = performance.now();
      if (st.at) showSpot(st.at, st.label);
      await runActs(st.do, run, false);
      if (st.after) showSpot(st.after, st.afterLabel);
      const used = (performance.now() - t0) * S.speed;
      await sleep(Math.max(st.do && st.do.length ? 1400 : 0, readMs(st.say) - used), run);
      hideSpot();
      check(run);
      if (S.si + 1 < c.steps.length) S.si++;
      else if (S.ci + 1 < CH.length) { jump(S.ci + 1, 0, false); return; }
      else { setPlaying(false); if (S.onEnd) S.onEnd(); return; }
    }
  }

  /* ---------------- controls ---------------- */
  function toast(m) { const t = el("div", "toast", esc(m)); document.body.append(t); setTimeout(() => t.remove(), 4200); }
  $("#bPlay").onclick = () => setPlaying(!S.playing);
  $("#bNext").onclick = () => step(1);
  $("#bPrev").onclick = () => step(-1);
  $("#bNextC").onclick = () => S.ci + 1 < CH.length && jump(S.ci + 1, 0, false);
  $("#bPrevC").onclick = () => jump(Math.max(0, S.si > 0 ? S.ci : S.ci - 1), 0, false);
  $("#bSpeed").onclick = () => { S.speed = SPEEDS[(SPEEDS.indexOf(S.speed) + 1) % SPEEDS.length]; $("#bSpeed").textContent = S.speed + "×"; };
  function step(d) {
    let ci = S.ci, si = S.si + d;
    if (si < 0) { if (ci === 0) return; ci--; si = CH[ci].steps.length - 1; }
    if (si >= CH[ci].steps.length) { if (ci + 1 >= CH.length) return; ci++; si = 0; }
    jump(ci, si, false);
  }
  const toggle = (btn, cls) => { $(btn).onclick = () => { $(btn).classList.toggle("on"); $("#app").classList.toggle(cls); fit(); }; };
  toggle("#bNav", "noNav"); toggle("#bCap", "noCap"); toggle("#bNotes", "noNotes");
  $("#bHelp").onclick = () => $("#help").showModal();
  document.addEventListener("keydown", (e) => {
    if (e.target.closest && e.target.closest("input,textarea,select")) return;
    if (e.code === "Space") { e.preventDefault(); setPlaying(!S.playing); }
    else if (e.key === "ArrowRight") { e.preventDefault(); e.shiftKey ? $("#bNextC").click() : step(1); }
    else if (e.key === "ArrowLeft") { e.preventDefault(); e.shiftKey ? $("#bPrevC").click() : step(-1); }
  });
  // the progress bar: click or drag to a step
  const track = $("#track"), tip = $("#hoverTip");
  const at = (ev) => { const r = track.getBoundingClientRect(); const f = Math.min(1, Math.max(0, (ev.clientX - r.left) / r.width)); const t = f * S.total; return { f, x: TL.find((s) => t < s.t0 + s.ms) || TL[TL.length - 1] }; };
  track.addEventListener("mousemove", (ev) => { const { f, x } = at(ev); tip.style.display = "block"; tip.style.left = f * 100 + "%"; tip.textContent = `${CH[x.ci].title} · ${CH[x.ci].steps[x.si].title}`; });
  track.addEventListener("mouseleave", () => { tip.style.display = "none"; });
  let drag = false;
  track.addEventListener("pointerdown", (ev) => { drag = true; track.setPointerCapture(ev.pointerId); const { f } = at(ev); $("#fill").style.width = f * 100 + "%"; $("#knob").style.left = f * 100 + "%"; });
  track.addEventListener("pointermove", (ev) => { if (!drag) return; const { f } = at(ev); $("#fill").style.width = f * 100 + "%"; $("#knob").style.left = f * 100 + "%"; });
  track.addEventListener("pointerup", (ev) => { if (!drag) return; drag = false; const { x } = at(ev); jump(x.ci, x.si, false); });

  /* ---------------- the stage scales to the space it has ---------------- */
  const APP_MIN_W = 1040;  // the app stacks its panels below 1000 px: it gets at least this width, scaled to fit
  function fitApp() {
    const win = $("#appwin"), sc = $("#appscale");
    const w = win.clientWidth, h = win.clientHeight, k = w < APP_MIN_W ? w / APP_MIN_W : 1;
    sc.style.width = w / k + "px"; sc.style.height = h / k + "px"; sc.style.transform = k === 1 ? "" : `scale(${k})`;
  }
  function fit() {
    fitApp();
    const wrap = $("#stageWrap"), stage = $("#stage");
    const sw = 1320, sh = 820, aw = wrap.clientWidth - 28, ah = wrap.clientHeight - 20;
    const k = Math.min(aw / sw, ah / sh);
    stage.style.transform = `scale(${k})`;
    stage.style.left = (wrap.clientWidth - sw * k) / 2 + "px";
    stage.style.top = Math.max(10, (wrap.clientHeight - sh * k) / 2) + "px";
  }
  window.addEventListener("resize", fit);

  /* ---------------- start ---------------- */
  (async () => {
    S.debug = new URLSearchParams(location.search).has("debug");
    S.embed = new URLSearchParams(location.search).has("embed");
    if (S.embed) { $("#app").classList.add("embed", "noNav", "noNotes"); }
    if (innerWidth < 1600) { $("#app").classList.add("noNav"); $("#bNav").classList.remove("on"); }
    const [m, a] = await Promise.all([fetch(DATA + "manifest.json").then((r) => r.json()), fetch(DATA + "agent.json").then((r) => r.json())]);
    S.data = m; S.agent = a;
    buildTimeline(); buildNav(); fit();
    $("#startSub").textContent = `${CH.length} chapters · ${Math.round(S.total / 60000)} min · pause anytime and try it yourself`;
    $("#start").onclick = () => setPlaying(true);
    const q = new URLSearchParams(location.search), ci = Math.max(0, CH.findIndex((c) => c.id === q.get("chapter")));
    $("#app").classList.add("paused");
    if (q.has("check")) {  // the developers' self-check: play every chapter fast, collect failed steps and unrecorded requests
      S.check = true; S.speed = Number(q.get("speed") || 3);
      const report = (what, data) => fetch("__check?" + new URLSearchParams({ what, data: JSON.stringify(data) })).catch(() => {});
      S.onEnd = () => { S.done = true; report("done", { fails: S.fails || [], misses: S.misses || [] }); };
      S.onStep = () => report("step", { ch: CH[S.ci].id, si: S.si, fails: (S.fails || []).length, misses: (S.misses || []).length });
      jump(ci, 0, true);
      return;
    }
    jump(ci, 0, false);
  })();
  window.__player = { S, jump, find, setScene, match, chapters: () => CH.map((c) => ({ id: c.id, title: c.title, ms: c.ms })) };  // for the developers' self-check (tools/tutorial/check.js)
})();
