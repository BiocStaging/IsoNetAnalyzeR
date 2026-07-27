with_temp_plot_device <- function(code) {
  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  device <- grDevices::dev.cur()
  on.exit({
    open_devices <- grDevices::dev.list()
    if (!is.null(open_devices) && device %in% open_devices) {
      grDevices::dev.off(device)
    }
    unlink(path)
  }, add = TRUE)
  force(code)
}
