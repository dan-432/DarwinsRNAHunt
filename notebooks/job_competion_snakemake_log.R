library(ggplot2)

log <- readLines("snakemake_6041397.log")

prog <- grep("^\\d+ of \\d+ steps \\(\\d+%\\) done$", log)
stamp <- grep("^\\[[A-Z][a-z]{2} ", log)

ts_chr <- rep(NA_character_, length(log))
ts_chr[stamp] <- sub("^\\[(.*)\\]$", "\\1", log[stamp])

last_ts <- NA_character_
for (i in seq_along(ts_chr)) {
  if (!is.na(ts_chr[i])) last_ts <- ts_chr[i]
  ts_chr[i] <- last_ts
}

df <- data.frame(
  time = as.POSIXct(ts_chr[prog], format = "%a %b %d %H:%M:%S %Y", tz = "UTC"),
  percent = as.integer(sub(".*\\((\\d+)%\\).*", "\\1", log[prog]))
)

df <- df[!is.na(df$time), ]

p <- ggplot(df, aes(x = time, y = percent)) +
  geom_line() +
  geom_point() +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) +
  labs(x = "Time", y = "Progress (%)", title = "Snakemake Progress Over Time") +
  theme_minimal()

ggsave("snakemake_progress.png", plot = p, width = 8, height = 4.5, dpi = 300)
