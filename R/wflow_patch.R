generate_rmd <- function(path, alias, dir) {
  path <- paste0(dir, '/', path)
  abs.path <- tools::file_path_as_absolute(paste0('../', path))
  cat(
    "---\n",
    yaml::as.yaml(rmarkdown::yaml_front_matter(abs.path)),
    "---\n\n",
    "**Source file\\:** ", path, "\n\n",
    "```{r child='", abs.path, "'}\n```",
    sep="",
    file=alias
  )
}

wflow_dir_build <- function(files=NULL, dir='notebooks', ...) {
  if (is.null(files)) {
    files <- list.files(paste0('../', dir), recursive=T, include.dirs=T, pattern="./*.(r|R)md")
  }
  else {
    for (file in files) {
      path <- paste0(dir, "/", file)
      if (!file.exists(paste0('../', path)))
        stop(paste0("File doesn't exist: ../", path))
    }
  }

  file_aliases <- gsub("/", "__", files)
  mapply(generate_rmd, files, file_aliases, dir=dir)
  workflowr::wflow_build(files=file_aliases, ...)
  invisible(file.remove(file_aliases))
}
