(function () {
  const cacheTtlMs = 60 * 60 * 1000;

  function compact(value) {
    const n = Number(value) || 0;
    if (n >= 1000) return `${(n / 1000).toFixed(n >= 10000 ? 0 : 1)}k`;
    return String(n);
  }

  function cacheKey(repo) {
    return `gw-repo-facts:${repo}`;
  }

  function readCache(repo) {
    try {
      const raw = sessionStorage.getItem(cacheKey(repo));
      if (!raw) return null;
      const parsed = JSON.parse(raw);
      if (!parsed || Date.now() - parsed.ts > cacheTtlMs) return null;
      return parsed.data;
    } catch {
      return null;
    }
  }

  function writeCache(repo, data) {
    try {
      sessionStorage.setItem(cacheKey(repo), JSON.stringify({ ts: Date.now(), data }));
    } catch {
      /* ignore quota */
    }
  }

  async function githubJson(url) {
    const response = await fetch(url, {
      headers: { Accept: "application/vnd.github+json" }
    });
    if (!response.ok) return null;
    return response.json();
  }

  async function loadFacts(repo) {
    const cached = readCache(repo);
    if (cached) return cached;

    const [meta, release] = await Promise.all([
      githubJson(`https://api.github.com/repos/${repo}`),
      githubJson(`https://api.github.com/repos/${repo}/releases/latest`)
    ]);
    if (!meta) return null;

    let tag = release && release.tag_name ? release.tag_name : "";
    if (!tag) {
      const tags = await githubJson(`https://api.github.com/repos/${repo}/tags?per_page=1`);
      if (Array.isArray(tags) && tags[0] && tags[0].name) tag = tags[0].name;
    }

    const data = {
      tag,
      stars: meta.stargazers_count,
      forks: meta.forks_count
    };
    writeCache(repo, data);
    return data;
  }

  function fillFact(root, suffix, text) {
    const el = root.querySelector(`.md-source__fact--${suffix}`);
    if (!el || !text) return;
    el.textContent = text;
    el.hidden = false;
  }

  async function hydrate(link) {
    const repo = link.getAttribute("data-gw-repo");
    if (!repo) return;
    const facts = await loadFacts(repo);
    if (!facts) return;
    fillFact(link, "version", facts.tag);
    fillFact(link, "stars", compact(facts.stars));
    fillFact(link, "forks", compact(facts.forks));
  }

  document.querySelectorAll(".gw-source[data-gw-repo]").forEach((link) => {
    hydrate(link);
  });
})();
