library(data.table)
library(dplyr)
library(ggpubr)
library(glmnet)
library(selectiveInference)
pair <- fread("03.EUR.tp.mr.info.res")
trna.info <- fread("01.lcl.eur.tRNA.classification.bed") %>%
  mutate(anticodon = str_split_i(ID,'-',3)) %>%
  group_by(anticodon) %>%
  reframe(N.isodecoder = n()) %>%
  mutate(cut_N.iso = cut(N.isodecoder,
                     breaks = c(0, 10, 20, 30),
                     labels = c("1-10", "11-20", "21-30"),
                     right = TRUE,
                     include.lowest = TRUE))


df.s = pair %>%
  mutate(anticodon = str_split_i(exposure,'-',3),
         b.adj = abs(b),
         quan.iv = ntile(as.numeric(nsnp),10)) %>%
  left_join(trna.info,by='anticodon')

ggplot(data= df.s, aes(y=abs(b),x = cut_N.iso,color=as.factor(cut_N.iso))) +
  geom_boxplot() +
  labs(x = 'Effect size') +
  scale_color_manual(values=c('#fcbba1','#ef3b2c','#99000d')) +
  # scale_color_manual(values=c('#99000d','#cb181d','#ef3b2c','#fc9272', '#fcbba1')) +
  theme_classic() +
  geom_signif(comparisons = list(c("1-10", "11-20"),c("1-10", "21-30")), 
              map_signif_level = FALSE,
              y_position = c(2.7, 3),
              color = 'black') +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))
ggsave("05.beta.description.plot(isodecoder.3group.boxplot).pdf",width = 5,height = 4)

