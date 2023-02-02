setwd(~/Github/DeCoTUR_manuscript_code)
#spydrpick runtime
# WITHOUT ARACNE
s100_1000 <- 2.375 + 0.524
s100_3000 <- 11.635+0.712
s100_5000 <- 29.719 + 1.005
s100_10000 <- 2*60+38.477+3
s500_1000 <- 5.543 + 0.525
s500_3000 <- 40.897 + 0.982
s500_5000 <- 1*60 + 51.583 + 1.205
s500_10000 <- 7*60+44.885+5.386
s1000_1000 <- 9.684+0.535
s1000_3000 <- 1*60+16.707+1.101
s1000_5000 <- 3*60+32.236 + 2.141
s1000_10000 <- 14*60+16.050+ 6.521
#

spydat <- c(s100_1000, s100_3000, s100_5000, s100_10000, 
            s500_1000, s500_3000, s500_5000, s500_10000,
            s1000_1000, s1000_3000, s1000_5000, s1000_10000)

spydat1 <- data.frame(sp = rep(c(100, 500, 1000), each = 4), 
                      gn = rep(c(1000, 3000, 5000, 10000), times = 3),
                      val = spydat)


ggplot(spydat1, aes(x = gn, y = val, group = sp, color = sp)) + geom_point() + 
  geom_line() + scale_color_viridis()
# Linear in sp, quad in gn

# WITH ARACNE
s100_1000a <- 22.775+0.601
s100_3000a <- 13*60+14.456+.870
s100_5000a <- 52*60+8.531+3.09
s100_10000a <- 810*60 + 47.947 + 22.181
s500_1000a <- 25.189 + 0.591
s500_3000a <- 16*60 + 37.868 + 1.09
s500_5000a <- 64*60 + 1.858 + 12.1
s500_10000a <- 725*60+44.055 + 32.420
s1000_1000a <- 39.878 + 0.571
s1000_3000a <- 20*60 + 33.711+2.612
s1000_5000a <- 86*60 + 17.860 + 18.011
s1000_10000a <- 638*60+44.132+31.612
#
spydat2 <- data.frame(sp=c(100, 100, 100, 100, 500, 500, 500, 500, 1000, 1000, 1000, 1000), 
                      gn = c(1000, 3000, 5000, 10000, 1000, 3000, 5000, 10000, 1000, 3000, 5000, 10000), 
                      val = c(s100_1000a, s100_3000a, s100_5000a, s100_10000a,
                              s500_1000a, s500_3000a, s500_5000a, s500_10000a,
                              s1000_1000a, s1000_3000a, s1000_5000a, s1000_10000a))


#spydat2 <- data.frame(sp = rep(c(100, 500, 1000), each = 4), 
#                      gn = rep(c(1000, 3000, 5000, 10000), times = 3),
#                      val = spydat)

#Coinfinder
c100_1000 <- 125.558 + 4.922 + 4*60+10.452 + 5.809
c100_3000 <- 12*60+1.333 + 14*60+9.699+6.165
c100_5000 <- 23*60+26.796+ 25*60+53.891+6.315
c100_10000 <- 26*60+25.035+6.649 + 23*60+51.759 + 6.411
c500_1000<- 32*60+20.312+12.178+33*60+13.390+12.285
c500_3000<- 85*60+25.313 + 54.345+77*60+21.664
c500_5000<- 112*60+5.718 + 100*60+55.285+21.781+25.838
c500_10000 <- 119*60 + 13.561 + 135*60+2.032+28.026
c1000_1000 <- 114*60+13.670+18.407+117*60+25.373+21.039
c1000_3000 <- 251*60+4.982+23.281+242*60+16.82+25.647
c1000_5000 <- 387*60+24.558+346*60+34.305+34.280
c1000_10000 <- 458*60+3.647 + 49.664 + 405*60+1.156 + 47.9
coindat <- c(c100_1000, c100_3000, c100_5000, c100_10000,
             c500_1000, c500_3000, c500_5000, c500_10000, 
             c1000_1000, c1000_3000, c1000_5000, c1000_10000)

coindat1 <-  data.frame(sp = rep(c(100, 500, 1000), each = 4), 
                       gn = rep(c(1000, 3000, 5000, 10000), times = 3),
                       val = coindat)

  
  
# Decotur 100
d100_1000_100 <- 7.86+0.89
d100_3000_100 <- 311.44
d100_5000_100 <- 1510.97
d100_10000_100 <- 13829.70
d500_1000_100 <- 19.86
d500_3000_100 <- 253.78
d500_5000_100 <- 1239.55
d500_10000_100 <- 15463.31  + 1691.28
d1000_1000_100 <- 60.33+1.25
d1000_3000_100 <- 351.41 + 22.16
d1000_5000_100 <- 1384.97 +  116.17
d1000_10000_100 <- 11115.61 + 2092.35

# Decotur 500
d100_1000_500 <- 43.31
d100_3000_500 <- 583.05
d100_5000_500 <- 2285.709
d100_10000_500 <- 17304.45
d500_1000_500 <- 50.55
d500_3000_500 <- 566.38
d500_5000_500 <- 2149.06
d500_10000_500 <- 16752.65
d1000_1000_500 <- 142.83 + 3.65
d1000_3000_500 <- 995.08 +  46.11
d1000_5000_500 <- 3419.75 + 174.32
d1000_10000_500 <- 16405.89 + 2406.5

# Decotur 1000
d100_1000_1000 <- 84.64
d100_3000_1000 <- 1055.16
d100_5000_1000 <- 3184.32
d100_10000_1000 <- 24253.23 +2916.7
d500_1000_1000 <- 103.03
d500_3000_1000 <- 950.90
d500_5000_1000 <- 3181.86
d500_10000_1000 <- 20436.42
d1000_1000_1000 <- 361.08 + 8.14
d1000_3000_1000 <- 2089.78  + 76.96
d1000_5000_1000 <-  7343.25 + 255.03
d1000_10000_1000 <- 27785.22+2959.3

# quadratic in cps and gn, constant in sp (as expected)

ncps <- c(rep(c(100, 500, 1000), each = 12))
gn <- rep(c(1000, 3000, 5000, 10000), times = 9)
sp <- rep(c(100, 500, 1000, 100, 500, 1000, 100, 500, 1000), each = 4)
decs <- c(d100_1000_100, d100_3000_100, d100_5000_100, d100_10000_100,
          d500_1000_100, d500_3000_100, d500_5000_100, d500_10000_100, 
          d1000_1000_100, d1000_3000_100, d1000_5000_100, d1000_10000_100,
          d100_1000_500, d100_3000_500, d100_5000_500, d100_10000_500, 
          d500_1000_500, d500_3000_500, d500_5000_500, d500_10000_500, 
          d1000_1000_500, d1000_3000_500, d1000_5000_500, d1000_10000_500,
          d100_1000_1000, d100_3000_1000, d100_5000_1000,  d100_10000_1000,
          d500_1000_1000, d500_3000_1000, d500_5000_1000, d500_10000_1000,
          d1000_1000_1000, d1000_3000_1000, d1000_5000_1000, d1000_10000_1000)
decdat <- data.frame(ncps, gn, sp, decs)
ggplot(decdat, aes(x = gn, y = decs, color = sp, group = sp)) + facet_wrap(~ncps) +
  scale_color_viridis() + geom_point() + geom_line()
  
spydat1$ncps <-NA
spydat1$type <- 'spy'
spydat2$type <- 'spya'
spydat2$ncps <- NA
decdat$type <- 'dec'
#head(spydat1)
decdat <- decdat[, c(3, 2, 4, 1, 5)]
names(decdat) <- c('sp', 'gn', 'val', 'ncps', 'type')
decdat$type <- paste0(decdat$ncps, 'dec')
coindat1$type <- 'coin'
coindat1$ncps <- NA
both <- rbind(spydat1, spydat2, decdat, coindat1)
mb <- melt(both, id = c('type', 'sp', 'gn', 'ncps'))
# timecomp <- ggplot(both, aes(x = gn, y = val, group = type, color = type)) + 
#   facet_wrap(~sp, labeller = labeller(sp = c('100' = '100 Samples', '500' = '500 Samples', '1000' = '1000 Samples'))) + 
#   geom_point() + 
#   geom_line(size=1) + 
#   scale_color_manual(name = 'Software', 
#                      labels = c(coin = 'Coinfinder', '100dec' =  'DeCoTUR (100 Close Pairs)', 
#                                 '1000dec' = 'DeCoTUR (1000 Close Pairs)', 
#                                 '500dec' = 'DeCoTUR (500 Close Pairs)', 
#                                 spy = 'SpydrPick (No ARACNE)', spya = 'SpydrPick (With ARACNE)'), 
#                      values = c('coin' = "#6DCD59FF", '1000dec' = "#440154FF",
#                                 '500dec' = "#414487FF",  '100dec' = "#2A788EFF",
#                                 'spya' ="#FDE725FF", 'spy' ="#B4DE2CFF")) + 
#   xlab('Number of Genes') + ylab('Runtime (CPU seconds)') + scale_y_log10() + theme_bw() + 
#   scale_x_log10(breaks = c(1000, 3000, 5000, 10000))
timecomp_loglin <- ggplot(both, aes(x = gn, y = val, group = type, color = type)) + 
  facet_wrap(~sp, labeller = labeller(sp = c('100' = '100 Samples', '500' = '500 Samples', '1000' = '1000 Samples'))) + 
  geom_point() + 
  geom_line(size=1) + 
  scale_color_manual(name = 'Software', 
                     labels = c(coin = 'Coinfinder', '100dec' =  'DeCoTUR (100 Close Pairs)', 
                                '1000dec' = 'DeCoTUR (1000 Close Pairs)', 
                                '500dec' = 'DeCoTUR (500 Close Pairs)', 
                                spy = 'SpydrPick (No ARACNE)', spya = 'SpydrPick (With ARACNE)'), 
                     values = c('coin' = "#6DCD59FF", '1000dec' = "#440154FF",
                                '500dec' = "#414487FF",  '100dec' = "#2A788EFF",
                                'spya' ="#FDE725FF", 'spy' ="#B4DE2CFF")) + 
  xlab('Number of Genes') + ylab('Runtime (CPU seconds)') + scale_y_log10() + theme_bw() + 
  scale_x_continuous(breaks = c(1000, 3000, 5000, 10000))
#save_plot('timecomp.pdf', timecomp, base_width = 12)
save_plot('figures/timecomp_loglin.pdf', timecomp_loglin, base_width = 12)
