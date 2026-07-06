/* Client-side filter + column sort for the Marker annotations table.
   Vanilla JS, no dependencies. Scoped to the marker-annotations page and
   re-run on Material's instant navigation via the document$ observable. */
document$.subscribe(function () {
  if (location.pathname.indexOf("marker-annotations") === -1) return;

  var table = document.querySelector(".md-content article table");
  if (!table || table.dataset.enhanced || !table.tBodies.length) return;
  table.dataset.enhanced = "1";

  var tbody = table.tBodies[0];
  var rows = Array.prototype.slice.call(tbody.rows);

  // Filter box above the table.
  var filter = document.createElement("input");
  filter.type = "search";
  filter.placeholder = "Filter " + rows.length + " markers by id or annotation…";
  filter.setAttribute("aria-label", "Filter markers");
  filter.style.cssText =
    "width:100%;box-sizing:border-box;padding:.5em .7em;margin:.4rem 0 .8rem;" +
    "border:1px solid var(--md-default-fg-color--lighter);border-radius:.2rem;" +
    "background:var(--md-default-bg-color);color:var(--md-default-fg-color);font-size:.75rem;";

  var anchor = table.closest(".md-typeset__scrollwrap") || table;
  anchor.parentNode.insertBefore(filter, anchor);

  filter.addEventListener("input", function () {
    var q = filter.value.trim().toLowerCase();
    for (var i = 0; i < rows.length; i++) {
      var hit = !q || rows[i].textContent.toLowerCase().indexOf(q) !== -1;
      rows[i].style.display = hit ? "" : "none";
    }
  });

  // Click a header to sort by that column (toggles ascending/descending).
  var headers = table.tHead ? table.tHead.rows[0].cells : [];
  Array.prototype.forEach.call(headers, function (th, idx) {
    th.style.cursor = "pointer";
    th.title = "Sort by " + th.textContent.trim();
    th.addEventListener("click", function () {
      var asc = th.dataset.asc !== "true";
      Array.prototype.forEach.call(headers, function (h) {
        delete h.dataset.asc;
      });
      th.dataset.asc = asc ? "true" : "false";
      rows.slice()
        .sort(function (a, b) {
          var x = a.cells[idx].textContent.trim().toLowerCase();
          var y = b.cells[idx].textContent.trim().toLowerCase();
          if (x < y) return asc ? -1 : 1;
          if (x > y) return asc ? 1 : -1;
          return 0;
        })
        .forEach(function (r) {
          tbody.appendChild(r);
        });
    });
  });
});
