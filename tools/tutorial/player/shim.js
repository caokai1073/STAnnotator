/* STAnnotator tutorial: replay shim. The build puts this first in the copied app page. Every request the app makes
 * to its server (/api/...) is answered from the recording of a real session, chosen by the tutorial player (the
 * parent page) for the current scene. The app itself is unchanged: the same page, the same code. */
(() => {
  "use strict";
  const T = window.parent && window.parent !== window ? window.parent.__TUT : null;
  if (!T) return;
  const isApi = (u) => u.origin === location.origin && (u.pathname.startsWith("/api/") || u.pathname.startsWith("/review/"));
  const NOT_IN_DEMO = "Not in the tutorial: this demo replays a recorded session. Install STAnnotator to run it on your own data.";

  // the app's own settings (last tab, opacity, drafts) live in memory: every chapter starts as a fresh app, and nothing
  // is written to the website's storage
  const mem = new Map();
  const storage = {
    getItem: (k) => (mem.has(String(k)) ? mem.get(String(k)) : null), setItem: (k, v) => { mem.set(String(k), String(v)); },
    removeItem: (k) => { mem.delete(String(k)); }, clear: () => mem.clear(), key: (i) => [...mem.keys()][i] ?? null,
    get length() { return mem.size; },
  };
  try { Object.defineProperty(window, "localStorage", { configurable: true, get: () => storage }); } catch (e) { /* keep the real one */ }

  // the browser's own confirm dialog (e.g. "Move the cells there onto their nuclei?") would stop the replay: the player
  // draws it for a moment, and the answer is yes, as in the recording
  window.confirm = (msg) => { T.confirm && T.confirm(String(msg || "")); return true; };

  const realFetch = window.fetch.bind(window);
  window.fetch = async (input, init) => {
    const url = new URL(typeof input === "string" ? input : input.url, location.href);
    if (!isApi(url)) return realFetch(input, init);
    const method = ((init && init.method) || (typeof input !== "string" && input.method) || "GET").toUpperCase();
    let body = init && init.body;
    if (body && typeof body !== "string") body = null;
    const hit = T.match(method, url.pathname + url.search, body || "");
    T.log && T.log(method, url.pathname + url.search, !!hit);
    if (!hit) {
      return new Response(JSON.stringify({ error: NOT_IN_DEMO }), { status: 409, headers: { "content-type": "application/json" } });
    }
    if (hit.delay) await new Promise((r) => setTimeout(r, hit.delay));
    const r = await realFetch(hit.href);
    return new Response(await r.arrayBuffer(), { status: hit.status, headers: { "content-type": hit.ctype, ...(hit.headers || {}) } });
  };

  // images (H&E tiles, figures): <img src="/api/..."> loads the recorded file instead
  const desc = Object.getOwnPropertyDescriptor(HTMLImageElement.prototype, "src");
  const mapSrc = (v) => {
    try {
      const u = new URL(String(v), location.href);
      if (!isApi(u)) return v;
      const hit = T.match("GET", u.pathname + u.search, "");
      return hit ? new URL(hit.href, location.href).href : T.blank;
    } catch (e) { return v; }
  };
  Object.defineProperty(HTMLImageElement.prototype, "src", {
    configurable: true, enumerable: desc.enumerable,
    get() { return desc.get.call(this); },
    set(v) { desc.set.call(this, mapSrc(v)); },
  });
  const setAttr = Element.prototype.setAttribute;
  Element.prototype.setAttribute = function (name, value) {
    if (this instanceof HTMLImageElement && String(name).toLowerCase() === "src") value = mapSrc(value);
    return setAttr.call(this, name, value);
  };
  // <img> written as HTML text (innerHTML) does not go through the setter: fix it when it enters the page
  const fixImgs = (root) => {
    const list = root.tagName === "IMG" ? [root] : root.querySelectorAll ? root.querySelectorAll("img[src^='/api/']") : [];
    for (const im of list) { const v = im.getAttribute("src"); if (v && v.startsWith("/api/")) desc.set.call(im, mapSrc(v)); }
  };
  new MutationObserver((ms) => { for (const m of ms) for (const n of m.addedNodes) if (n.nodeType === 1) fixImgs(n); })
    .observe(document.documentElement, { childList: true, subtree: true });
  // downloads and new windows stay inside the tutorial
  window.open = () => { T.toast && T.toast(NOT_IN_DEMO); return null; };
  document.addEventListener("click", (e) => {
    const a = e.target.closest && e.target.closest("a[href]");
    if (!a) return;
    const u = new URL(a.getAttribute("href"), location.href);
    if (isApi(u) || a.target === "_blank") { e.preventDefault(); T.download ? T.download(u.pathname + u.search, a) : null; }
  }, true);
  T.attach && T.attach(window);
})();
