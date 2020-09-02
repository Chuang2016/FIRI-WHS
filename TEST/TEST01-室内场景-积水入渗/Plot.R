library(data.table)
library(ggplot2)

data <- fread(file = "SWF_PER.OUT", header = T, skip = 11)
ggplot(data, aes(x = TIME, y = WTOP/PERLEN)) +
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
ggplot(data, aes(x = MOIST, y = DEPTH, color = factor(TIME))) +
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

ggplot(data, aes(x = TIME, y = TEMP, color = factor(DEPTH))) +
  theme_bw() + 
  geom_line() +
  xlab("时间 (min)") +
  ylab("温度 (°C)")

data <- fread(file = "OBS.OUT", header = T, skip = 9)
ggplot(data, aes(x = TIME, y = CONC, color = factor(DEPTH))) +
  theme_bw() + 
  geom_line() +
  xlab("时间 (min)") +
  ylab("浓度 (mg/cm3)")
.
