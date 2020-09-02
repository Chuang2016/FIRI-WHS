library(data.table)
library(ggplot2)

data <- fread(file = "SWF_PER.OUT", header = T, skip = 11)
ggplot(data, aes(x = TIME, y = WBOT/PERLEN)) +
  theme_bw() + 
  geom_line() +
  xlab("时间 (min)") + 
  ylab("地表水分入渗速率 (cm/min)")

data <- fread(file = "SWF_PER.OUT", header = T, skip = 11)
ggplot(data, aes(x = TIME, y = C_WTOP)) +
  theme_bw() + 
  geom_line() +
  xlab("时间 (min)") + 
  ylab("地表水分入渗累计量 (cm)")


data <- fread(file = "SWF_PFL.OUT", header = T, skip = 12)
ggplot(data[TIME==120], aes(x = MOIST, y = DEPTH, color = factor(TIME))) +
  theme_bw() + 
  geom_point() +
  scale_y_reverse()

data <- fread(file = "SWF_OBS.OUT", header = T, skip = 11)
ggplot(data, aes(x = TIME, y = HEAD, color = factor(DEPTH))) +
  theme_bw() + 
  geom_line() +
  xlab("时间 (min)") +
  ylab("水势 (cm)")

ggplot(data, aes(x = TIME, y = MOIST, color = factor(DEPTH))) +
  theme_bw() + 
  geom_line() +
  xlab("时间 (min)") +
  ylab("含水率 (cm3/cm3)")
