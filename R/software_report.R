# Software citations report ----

# Add navigation to the rendered software report by linking to the headings
# already present in the HTML.
add_software_report_contents <- function(path) {
  html <- readLines(path, encoding = "UTF-8", warn = FALSE)
  if (any(grepl('id="software-contents"', html, fixed = TRUE))) return(path)
  package_section <- grep('<div id="r-packages-used"', html, fixed = TRUE)
  citations_section <- grep('<div id="package-citations"', html, fixed = TRUE)
  if (length(package_section) != 1L || length(citations_section) != 1L ||
      package_section >= citations_section) {
    stop("Could not find the software report sections for the contents links.")
  }

  contents <- c(
    '<nav id="software-contents" aria-label="Contents">',
    '<h2>Contents</h2>',
    '<ul>',
    '<li><a href="#r-packages-used">R packages used</a></li>',
    '<li><a href="#package-citations">Package citations</a></li>',
    '</ul>',
    '</nav>'
  )

  writeLines(enc2utf8(append(html, contents, after = package_section - 1L)),
             path, useBytes = TRUE)
  path
}
