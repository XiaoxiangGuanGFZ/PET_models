library(shiny)
library(ggplot2)
library(ggthemes)
library(shinythemes)
library(RColorBrewer)
library(Cairo)   # For nicer ggplot2 output when deployed on Linux
#----------------全局函数------------

Penman = function(Ta,Tmax,Tmin,RH,Wv,Pa,n,J,lat) { 
  #daily model
  
  #依据日平均气温Ta（℃），相对湿度RH（%），风速Wv（m3/s）
  #日照时数h(小时),大气压Pa(kPa,千帕),基于彭曼蒸散发公式
  #计算日潜在蒸散发，mm。
  
  # J:时序，一年的第多少天
  # lat:纬度，单位为：°，在计算中需要转换成弧度制（rad）
  # n:实际日照时数，小时
  
  e0 = function(Ta) {
    #根据气温，计算得到瞬时饱和水汽压kPa
    return(0.6108*exp(17.277*Ta/(Ta+237.3)))
  }
  
  u2 = Wv*4.87/(log(67.8*10-5.42))  #u2表示2m高的风速；
  gama = 0.665*10^(-3)*Pa     #gama为干湿常数，单位：kPa*℃^(-1)
  es = 0.5*(e0(Tmax)+e0(Tmin)) #es为饱和水汽压
  ea = RH*es/100  #ea为实际水汽压，kPa
  delta = 4098*e0(Ta)/(Ta+237.3)^2  #delta为饱和水汽压曲线斜率（kPa*℃^-1）
  
  #****计算地球外辐射，MJ/m2/d********
  dr =1+0.033*cos(2*pi/365*J)  #反转日地平均距离
  del = 0.408*sin(2*pi/365*J-1.39)  #太阳磁偏角（rad）
  w_s = acos(-tan(lat*pi/180)*tan(del))    #日出时角(rad)
  Ra = 37.59*dr*(w_s*sin(lat*pi/180)*sin(del)+cos(lat*pi/180)*cos(del)*sin(w_s)) #Ra：地球外辐射
  
  N = 24/pi*w_s #N：最大可能日照时数，小时
  Rs0 = (0.25+0.5)*Ra
  Rs = (0.25+0.5*n/N)*Ra
  
  Rns = (1-0.23)*Rs
  Rnl = 4.903*10^(-9)*((Tmax+273.16)^4+(Tmin+273.16)^4)*0.5*(0.34-0.14*sqrt(ea))*(1.35*Rs/Rs0-0.35)
  Rn = Rns-Rnl  #Rn:净辐射，MJ/m2/d
  
  #********************************************
  ET = (0.408*delta*Rn+gama*900/(Ta+273)*u2*(es-ea))/(delta+gama*(1+0.34*u2))
  return(round(ET,digits = 2))
  
}


Hamon = function(D,Ta,ha) { 
  #daily model
  #ha:常数，0.55
  #D 为日照时数,输入时单位为h，在公式里换算成12h；
  #Ta 为日平均气温，℃
  ET= ha*25.4*(D/12)^2*4.95*exp(0.062*Ta)/100
  return(round(ET,digits = 2))
}

Thornthwaite = function(C,I,Ta,d,N) {
  #daily model
  
  #C:为公式中一个常数，可以依据实际计算站点而变化，在公式的提出区域取值为16
  #I：年热力指数，提前算好，导入到函数中
  #d：日的日照时间，以小时计；
  #N：所在月的总天数，比如一月，N为31；6月份，N为30
  a = 67.5*10^(-8)*I^3-77.1*10^(-6)*I^2+0.0179*I+0.492
  ET = C*(10*Ta/I)^a*(d/12)*N/30
  return(round(ET,digits = 2))
  
}

Linacre = function(Ta,lat,RH,L) {
  #daily model
  #Ta:日平均气温，℃
  #Td:露点温度，℃
  #L :a constant, 500
  Td = 243.04*(log(RH/100)+17.625*Ta/(243.04+Ta))/(17.625-log(RH/100)-17.625*Ta/(243.04+Ta))
  ET= (L*(Ta+0.006*2)/(100-lat)+15*(Ta-Td))/(80-Ta)
  return(round(ET,digits = 2))
}


Blaney_Criddle = function(k,Ta,p) {
  #daily method
  
  #k:月消耗系数，可能每个月不一样，也可以全年用一个值，通常取值范围为[0.5,1.2]，常取0.85;
  #k值可以依据实际站点而变化；
  
  #Ta:日平均气温，℃
  #p:一天当中的日照时数与全年日照时数的比例，无量纲，范围为[0,1]
  ET = k*p*(0.46*Ta+8.13)
  return(round(ET,digits = 2))
  
}

Kharrufa = function(Ta,p,kh) {
  #daily method
  #kh: a constant,0.34
  #Ta:日平均气温，℃
  #p:一天当中的日照时数与全年日照时数的比例，无量纲，范围为[0,1]
  
  ET = kh*p*Ta^1.3
  return(round(ET,digits = 2))
}

Hargreaves = function(k,Ta,RH,lat,J,Tmax,Tmin) {
  #daily model
  # J:时序，一年的第多少天
  # lat:纬度，单位为：°，在计算中需要转换成弧度制（rad）
  
  # k:为一个常系数，可以变动以适应计算的站点,常取值为0.0023
  #****计算地球外辐射，MJ/m2/d********
  dr =1+0.033*cos(2*pi/365*J)  #反转日地平均距离
  del = 0.408*sin(2*pi/365*J-1.39)  #太阳磁偏角（rad）
  w_s = acos(-tan(lat*pi/180)*tan(del))    #日出时角(rad)
  Ra = 37.59*dr*(w_s*sin(lat*pi/180)*sin(del)+cos(lat*pi/180)*cos(del)*sin(w_s)) #Ra：地球外辐射
  
  TD = Tmax-Tmin
  ET = k*Ra*TD^0.5*(Ta+17.8)
  return(round(ET,digits = 2))
  
}

FAO24 = function(Ta,Pa,n,J,lat,c) { #Doorenbos-Pruitt (1977)
  #daily model
  # c:调整系数，与平均湿度，日照时间，风速有关。
  #依据日平均气温Ta（℃）
  #日照时数h(小时),大气压Pa(kPa,千帕),基于彭曼蒸散发公式
  #计算日潜在蒸散发，mm。
  
  # lat:纬度，单位为：°，在计算中需要转换成弧度制（rad）
  # n:实际日照时数，小时
  
  e0 = function(Ta) {
    #根据气温，计算得到瞬时饱和水汽压kPa
    return(0.6108*exp(17.277*Ta/(Ta+237.3)))
  }

  gama = 0.665*10^(-3)*Pa     #gama为干湿常数，单位：kPa*℃^(-1)

  delta = 4098*e0(Ta)/(Ta+237.3)^2  #delta为饱和水汽压曲线斜率（kPa*℃^-1）
  
  #****计算地球外辐射，MJ/m2/d********
  dr =1+0.033*cos(2*pi/365*J)  #反转日地平均距离
  del = 0.408*sin(2*pi/365*J-1.39)  #太阳磁偏角（rad）
  w_s = acos(-tan(lat*pi/180)*tan(del))    #日出时角(rad)
  Ra = 37.59*dr*(w_s*sin(lat*pi/180)*sin(del)+cos(lat*pi/180)*cos(del)*sin(w_s)) #Ra：地球外辐射
  
  N = 24/pi*w_s #N：最大可能日照时数，小时
  Rs0 = (0.25+0.5)*Ra
  Rs = (0.25+0.5*n/N)*Ra
  
  
  #********************************************
  ET = c*delta/(delta+gama)*Rs
  return(round(ET,digits = 2))
  
}


Makkink = function(Ta,Pa,n,J,lat) { #1957,无参数
  #daily model

  #依据日平均气温Ta（℃）
  #日照时数h(小时),大气压Pa(kPa,千帕),基于彭曼蒸散发公式
  #计算日潜在蒸散发，mm。
  
  # lat:纬度，单位为：°，在计算中需要转换成弧度制（rad）
  # n:实际日照时数，小时
  
  e0 = function(Ta) {
    #根据气温，计算得到瞬时饱和水汽压kPa
    return(0.6108*exp(17.277*Ta/(Ta+237.3)))
  }
  
  gama = 0.665*10^(-3)*Pa     #gama为干湿常数，单位：kPa*℃^(-1)
  
  delta = 4098*e0(Ta)/(Ta+237.3)^2  #delta为饱和水汽压曲线斜率（kPa*℃^-1）
  
  #****计算地球外辐射，MJ/m2/d********
  dr =1+0.033*cos(2*pi/365*J)  #反转日地平均距离
  del = 0.408*sin(2*pi/365*J-1.39)  #太阳磁偏角（rad）
  w_s = acos(-tan(lat*pi/180)*tan(del))    #日出时角(rad)
  Ra = 37.59*dr*(w_s*sin(lat*pi/180)*sin(del)+cos(lat*pi/180)*cos(del)*sin(w_s)) #Ra：地球外辐射
  
  N = 24/pi*w_s #N：最大可能日照时数，小时
  Rs0 = (0.25+0.5)*Ra
  Rs = (0.25+0.5*n/N)*Ra
  
  
  #********************************************
  ET = 0.61*delta/(delta+gama)*Rs/2.45-0.12
  return(round(ET,digits = 2))
  
}


Romanenko = function(Ta,RH,ro) {
  #monthly method
  #Ta,月平均气温，℃；
  #RH，月相对湿度，1%
  #ro: a constant 0.0018
  ET = ro*(25+Ta)^2*(100-RH)
  return(round(ET,digits = 2))
  
}

#---------Turc method---------
Turc = function(Ta,RH,Pa,n,J,lat) {
  r = 2.45  #汽化潜热MJ/kg
  At = NULL
  for (i in 1:length(RH)) {
    if (RH[i] >= 50) {
      At[i] = 1
    } else {
      At[i] = 1+(50-RH[i])/70
    }
  }
  #****计算太阳辐射，MJ/m2/d********
  dr =1+0.033*cos(2*pi/365*J)  #反转日地平均距离
  del = 0.408*sin(2*pi/365*J-1.39)  #太阳磁偏角（rad）
  w_s = acos(-tan(lat*pi/180)*tan(del))    #日出时角(rad)
  Ra = 37.59*dr*(w_s*sin(lat*pi/180)*sin(del)+cos(lat*pi/180)*cos(del)*sin(w_s)) #Ra：地球外辐射
  N = 24/pi*w_s #N：最大可能日照时数，小时
  Rs0 = (0.25+0.5)*Ra
  Rs = (0.25+0.5*n/N)*Ra
  
  #**********************
  PE = At*0.013*Ta/(Ta+15)*(23.885*Rs+50)/r
  return(round(PE,digits = 2))
}

#---------------Priestley-Taylor method-----------
Priestley_Taylor = function(Ta,Tmax,Tmin,RH,Pa,n,J,lat,a){
  #n:实际日照时数；
  #Ta：日平均气温；
  #Pa：气压，KPa
  #J：日序
  #lat：纬度
  #a:经验常数：1.26
  r = 2.45  #汽化潜热MJ/kg
  delta = 4098*e0(Ta)/(Ta+237.3)^2  #delta为饱和水汽压曲线斜率（kPa*℃^-1）
  
  e0 = function(Ta) {
    #根据气温，计算得到瞬时饱和水汽压kPa
    return(0.6108*exp(17.277*Ta/(Ta+237.3)))
  }
  
  gama = 0.665*10^(-3)*Pa     #gama为干湿常数，单位：kPa*℃^(-1)
  es = 0.5*(e0(Tmax)+e0(Tmin)) #es为饱和水汽压
  ea = RH*es/100  #ea为实际水汽压，kPa
  #****计算地球外辐射，MJ/m2/d********
  dr =1+0.033*cos(2*pi/365*J)  #反转日地平均距离
  del = 0.408*sin(2*pi/365*J-1.39)  #太阳磁偏角（rad）
  w_s = acos(-tan(lat*pi/180)*tan(del))    #日出时角(rad)
  Ra = 37.59*dr*(w_s*sin(lat*pi/180)*sin(del)+cos(lat*pi/180)*cos(del)*sin(w_s)) #Ra：地球外辐射
  N = 24/pi*w_s #N：最大可能日照时数，小时
  Rs0 = (0.25+0.5)*Ra
  Rs = (0.25+0.5*n/N)*Ra
  
  Rns = (1-0.23)*Rs
  Rnl = 4.903*10^(-9)*((Tmax+273.16)^4+(Tmin+273.16)^4)*0.5*(0.34-0.14*sqrt(ea))*(1.35*Rs/Rs0-0.35)
  Rn = Rns-Rnl  #Rn:净辐射，MJ/m2/d
  #**********************
  PE = a*delta/(delta+gama)*Rn/r
  return(round(PE,digits = 2))
  
}

DayofID = function(y,m,d) {
  run = c(31,29,31,30,31,30,31,31,30,31,30,31)
  ping = c(31,28,31,30,31,30,31,31,30,31,30,31)
  if ( (y%%4 == 0) & (y %%100 != 0) |(y %% 400 == 0) ) {   
    days=run
  } else {days=ping}
  
  if (m == 1) {
    J = d
  } else {
    J = sum(days[1:(m-1)])+d
  }
  return(J)
}


ui <- tagList(
  shinythemes::themeSelector(),
  navbarPage("PE",
             
             # ---------------Introduction------------------
             tabPanel("Potential Evapotranspiration",
                      fluidRow(column(12,
                                      wellPanel(
                                        h3("SOFTWARE PREFACE",align = "center")       
                                      ))),
                      fluidRow(column(12,
                                      wellPanel(
                                        p("水资源是人类生存和发展不可或缺的自然资源，目前水资源短缺已成为全世界关注的热点问题，而蒸散发既是水量平衡和地表热量的组成部分，又是影响水循环的重要因子。因此，蒸散发对水资源的形成过程、变化规律和水资源评价方面都有重要作用。用于表述蒸发变量的主要有水面蒸发、潜在蒸发和实际蒸散发。而潜在蒸散发（PE）是参考作物蒸散发即充分供水情况下自由水面或下垫面有植被覆盖的区域的蒸散发能力，能够全面地反映一个地区的蒸散发能力，是实际蒸散发的上限值。潜在蒸散发（PE）在作物需水量、水资源合理配置、气候干湿状况分析以及灌溉决策制定等方面研究得到广泛得重视与应用。
本软件提供了12种不同的潜在蒸散发模型，每种模型的计算原理、数据要求和使用情况各不相同，为广大科研工作者提供了方便的计算工具。
                                          使用过程中，如有问题请与软件开发者邮件联系，作者保留所有权利。
                                          "),
                                        br(),
                                        p("Programmer:Shane Guan",align = "right"),
                                        p("E-mail:xxguan@hhu.edu.cn",align = "right")  
                                        
                                        ))),
                      p(em("Address:College of hydrology and water resources, Hohai University,No.1,Xikang Road,NanJing 210098, China"),align = "center")
                      ),
             #--------------数据准备-----------------------
             tabPanel("模型选择",
                     "Description for different potential evapotranspiration methods",
                               navlistPanel(
                                 "Methodology",
                                 #-----------------
                                 tabPanel("Penman-Monteith method",
                                          fluidPage(
                                            h3("Daily method"),
                                            p("世界粮农组织（FAO）推荐的Penman-Monteith潜在蒸散发模型以水汽扩散理论和能量平衡为基础，
                                              相对全面的考虑了影响潜在蒸散发的各个因素，是目前世界上工人的在干旱区和湿润区精度较高，相对误差较小的方法。"),
                                            br(),
                                            p("Input elements:"),
                                            p("日平均、最高、最低气温(℃);日照时数 (小时);气压 (kPa,千帕);相对湿度(%);风速(m3/s);站点纬度(°)"),                        
                                            br(),
                                            p("For citation and detailed reference:"),
                                            p("Hamon W R. Estimating Potential Evapotranspiration[J]. Journal of Hydraulic Engineering, 1960, 87(3): 107-120.")
                                            
                                          )
                                          
                                 ),
                                 #---------------------------------------
                                 tabPanel("Hamon model",
                                          fluidPage(
                                            h3("Daily method"),
                                            p("Hamon公式可以计算得到各站点的逐日潜在蒸散发量。该方法的优点在于只要求温度数据，就能获得近似于复杂方法得到的结果。"),
                                            br(),
                                            p("Input elements:"),
                                            p("日平均气温(℃);日照时数 (小时);"),                        
                                            br(),
                                            p("For citation and detailed reference:"),
                                            p("R. D. Allen;L. A. S. Pereira;D. Raes;M. Smith.Crop evapotranspiration : guidelines for computing crop water requirements.FAO.1980,300,D05109")
                                            
                                          )
                                 ),
                                
                                 tabPanel("Thornthwaite method",
                                          fluidPage(
                                            h3("Daily method"),
                                            p("Thornthwaite (1948)基于水量平衡于水汽充足的山谷地区建立月平均气温与蒸散发的关系.需要说明的是，Thornthwaite潜在蒸散发模型是经验性再分配式蒸散发模型，在再分配过程考虑日尺度或者月尺度，这可能会对计算精度产生影响。"),
                                            br(),
                                            p("Input elements:"),
                                            p("日平均气温(℃);日照时数 (小时);"),                        
                                            br(),
                                            p("For citation and detailed reference:"),
                                            p("Thornthwaite C W. An Approach toward a Rational Classification of Climate[J]. Geographical Review, 1948, 38(1): 55-94")
                                            
                                          )
                                 ),
                                 tabPanel("Linacre method ",
                                          fluidPage(
                                            h3("Daily method"),
                                            p("在水分充足的植被环境，反射率约为0.25的情形下，Linacre（1977）[6]简化了Penman公式。"),
                                            br(),
                                            p("Input elements:"),
                                            p("日平均气温(℃);纬度(°);相对湿度(%)"),                        
                                            br(),
                                            p("For citation and detailed reference:"),
                                            p("Linacre E T. A simple formula for estimating evaporation rates in various climates, using temperature data alone[J]. Agricultural Meteorology, 1977, 18(6): 409-424.")
                                            
                                          )
                                 ),
                                 tabPanel("Blaney-Criddle method ",
                                          fluidPage(
                                            h3("Daily method"),
                                            p("美国西部的经验性公式。"),
                                            br(),
                                            p("Input elements:"),
                                            p("日平均气温(℃);日照时数 (小时);"),                        
                                            br(),
                                            p("Note:当温度低于-8.13℃时，PE计算值为0"),
                                            p("For citation and detailed reference:"),
                                            p(" Blaney H F, Criddle W D. Determining Water Requirements in Irrigated Areas From Climatological and Irrigation Data[J], 1952.")
                                            
                                          )
                                 ),
                                 tabPanel("Kharrufa method ",
                                          fluidPage(
                                            h3("Daily method"),
                                            p("Kharrufa(1985)依据PE/P与T的关系，得到公式。"),
                                            br(),
                                            p("Input elements:"),
                                            p("日平均气温(℃);日照时数 (小时);"),                        
                                            br(),
                                            p("Note:当温度低于0℃时，PE计算值为0"),
                                            p("For citation and detailed reference:"),
                                            p(" Xu C Y, Singh V P. Evaluation and generalization of temperature-based methods for calculating evaporation[J]. Hydrological Processes, 2001, 15(2): 339-349.")
                                            
                                          )
                                 ),
                                 tabPanel("Hargreaves method ",
                                          fluidPage(
                                            h3("Daily method"),
                                            p("由于太阳辐射数据通常情况下难以获得，Hargreaves (1985)用地球外辐射代替太阳辐射计算潜在蒸散发量"),
                                            br(),
                                            p("Input elements:"),
                                            p("日平均、最高、最低气温(℃);纬度(°);相对湿度(%)"),                        
                                            br(),
                                            p("Note:当平均温度低于-17.8℃时，PE计算值为0"),
                                            p("For citation and detailed reference:"),
                                            p(" H. Hargreaves G, Samani Z. Reference Crop Evapotranspiration From Temperature[M]. 1.1985.")
                                            
                                          )
                                 ),
                                 tabPanel(" Romanenko method ",
                                          fluidPage(
                                            h3("Monthly method"),
                                            p("Romanenko (1961)依据平均气温与相对湿度的关系建立潜在蒸散发模型。"),
                                            br(),
                                            p("Input elements:"),
                                            p("monthly平均气温(℃);monthly相对湿度(%)"),                        
                                            br(),
                                            p("For citation and detailed reference:"),
                                            p("Romanenkov.A. Computation of the autumn soil moisture using a universal relationship for a large area[J], 1961.")
                                            
                                          )
                                 ),
                                 tabPanel("FAO-24 Radiation method",
                                          fluidPage(
                                            h3("Daily method"),
                                            p("基于辐射的模型:Doorenbos and Pruitt (1977):其中太阳辐射依据理论公式计算得出，并非实测值。"),
                                            br(),
                                            p("Input elements:"),
                                            p("daily平均气温(℃);日照时数(h);气压(Pa)"),                        
                                            br(),
                                            p("For citation and detailed reference:"),
                                            p("")
                                            
                                          )
                                 ),
                                 tabPanel("Makkink method",
                                          fluidPage(
                                            h3("Daily method"),
                                            p("基于辐射的模型:Makkink (1957):其中太阳辐射依据理论公式计算得出，并非实测值。"),
                                            br(),
                                            p("Input elements:"),
                                            p("daily最高、最低气温(℃);日照时数(h);气压(Pa)"),                        
                                            br(),
                                            p("For citation and detailed reference:"),
                                            p("")
                                            
                                          )
                                 ),
                                 tabPanel("Priestley-Taylor method",
                                          fluidPage(
                                            h3("Daily method"),
                                            p("基于辐射的模型:Priestley C.H.B, Taylor R.J. (1972):其中太阳辐射依据理论公式计算得出，并非实测值。"),
                                            br(),
                                            p("Input elements:"),
                                            #***Ta,Tmax,Tmin,RH,Pa,n,J,lat,
                                            p("daily平均、最高、最低气温(℃);相对湿度（%）;日照时数(h);气压(Pa)"),                        
                                            br(),
                                            p("For citation and detailed reference:"),
                                            p("Priestley C H B, Taylor R J. On the Assessment of Surface Heat Flux and Evaporation Using Large-Scale Parameters[J]. Monthly Weather Review, 1972, 100(2): 81-92.")
                                            
                                          )
                                 ),
                                 tabPanel("Turc method",
                                          fluidPage(
                                            h3("Daily method"),
                                            p("基于辐射的模型:Turc(1961):其中太阳辐射依据理论公式计算得出，并非实测值。"),
                                            br(),
                                            p("Input elements:"),
                                            #***Ta,Tmax,Tmin,RH,Pa,n,J,lat,
                                            p("daily最高、最低气温(℃);相对湿度（%）;日照时数(h);气压(Pa)"),                        
                                            br(),
                                            p("For citation and detailed reference:"),
                                            p("")
                                            
                                          )
                                 )
                                 
                                 
                                 
                               ),
                     fluidRow(
                       column(6,
                              wellPanel(
                                #-----PE模型选择-------
                                selectInput("PEmodel", 
                                            "Please choose a PE model for calculation:",
                                            choices = c("Penman-Monteith","Hamon","Thornthwaite","Linacre",
                                                        "Blaney_Criddle","Kharrufa","Hargreaves","Romanenko",
                                                        "FAO-24 Radiation","Makkink","Priestley-Taylor","Turc")
                                            ,selected = "Penman-Monteith")
                                
                              )),
                       column(6,wellPanel(
                         conditionalPanel("input.PEmodel == 'Penman-Monteith' | input.PEmodel == 'Linacre' | input.PEmodel == 'Hargreaves' | input.PEmodel == 'FAO-24 Radiation' | input.PEmodel == 'Makkink' | input.PEmodel == 'Priestley-Taylor' | input.PEmodel == 'Turc' ",
                                       numericInput("lat","请输入站点纬度(北纬为正，南纬为负；单位：°):",value = 24.5,step = 0.01)        
                          ),
                         p("Please go to input data!")
                       ))
                     )
                     
                      
             ),
             #-----------资料输入模块------------------
             tabPanel("资料输入",
                      fluidRow(
                        column(4,
                               wellPanel(
                                 
                                 radioButtons("Sequencetype", "Sequence type",
                                              choices = c(
                                                monthly = "monthly",
                                                daily = "daily"),
                                             selected = "daily"),
                                 selectInput("Null", 
                                             "Null value:",
                                             choices = c(-99.9,-99,NA)
                                             ,selected = -99),
                                 tags$hr(),
                                 checkboxInput("header", "Header", TRUE),
                                 radioButtons("sep", "Separator",
                                              choices = c(Comma = ",",
                                                          Semicolon = ";",
                                                          Tab = "\t"),
                                              selected = ","),
                                 radioButtons("quote", "Quote",
                                              choices = c(None = "",
                                                          "Double Quote" = '"',
                                                          "Single Quote" = "'"),
                                              selected = '"')
                                 
                               )),
                        column(8,wellPanel(
                          fluidRow(
                            #-------导入文件-----
                            column(12,wellPanel(
                              fileInput("DATA", "Choose data File",
                                        multiple = TRUE,
                                        accept = c("text/csv",
                                                   "text/comma-separated-values,text/plain",
                                                   ".csv"))
                            ))
                          ),
                          fluidRow(h4("Data sequence range:",align = "center")),
                          fluidRow(
                            column(12, wellPanel(tableOutput("range")))
                        
                          ),
                          fluidRow(h4("Data sequence Display:",align = "center")),
                          fluidRow(
                            column(12,wellPanel(
                              DT::dataTableOutput("DATA_show")
                            ))
                          )
                        ))
                      )
                      
             ),
             
             #----------潜在蒸散发计算------------------
             tabPanel("蒸散发计算",
                      fluidRow(
                        column(4,h4("计算结果"),
                               #----------条件面板：模型参数--------
                               sliderInput("period", "Select baseline period", min = 1950, 
                                           max = 2025, value = c(1965, 1980),step = 1),
                               conditionalPanel("input.PEmodel == 'Thornthwaite'",
                                                numericInput("C","Constant C:",value = 16,min = 2,step = 2)),
                               conditionalPanel("input.PEmodel == 'Linacre'",
                                                numericInput("L","Constant L:",value = 500,min = 2,step = 2)),
                               conditionalPanel("input.PEmodel == 'Hargreaves'",
                                                numericInput("Ha","Constant H:",value = 0.0023,min = 0,step = 0.0002)),
                               conditionalPanel("input.PEmodel == 'Blaney_Criddle'",
                                                numericInput("k","Constant k:",value = 0.85,min = 0,step = 0.1)),
                               conditionalPanel("input.PEmodel == 'Kharrufa'",
                                                numericInput("kh","Constant kh:",value = 0.34,min = 0.1,step = 0.02)),
                               conditionalPanel("input.PEmodel == 'Hamon'",
                                                numericInput("ha","Constant ha:",value = 0.55,min = 0.01,step = 0.02)),
                               conditionalPanel("input.PEmodel == 'Romanenko'",
                                                numericInput("ro","Constant ro:",value = 0.0018,min = 0,step = 0.00001)),
                               conditionalPanel("input.PEmodel == 'FAO-24 Radiation'",
                                                numericInput("c","Constant c:",value = 1,min = 0,step = 0.005)),
                               conditionalPanel("input.PEmodel == 'Priestley-Taylor'",
                                                numericInput("alfa","Constant alfa:",value = 1.26,min = 0,step = 0.05)),
                               br(),
                               br(),
                               br(),
                               radioButtons("Real_PE", "是否有蒸发资料作比较？",
                                            choices = c(Yes = "Yes",
                                                        No = "No"),
                                            selected = "No"),
                               conditionalPanel("input.Real_PE == 'Yes' ",
                                                radioButtons("Sequencetype2", "Sequence type",
                                                             choices = c(
                                                               monthly = "monthly",
                                                               daily = "daily"),
                                                             selected = "daily"),
                                                fileInput("RealPE", "Choose the reference data",
                                                          multiple = TRUE,
                                                          accept = c("text/csv",
                                                                     "text/comma-separated-values,text/plain",
                                                                     ".csv"))
                                                
                               ),
                               #----计算结果下载-----
                               downloadButton("downloadData","点击下载潜在蒸散发计算结果")
                               
                        ),
                        column(8,
                               #----------PE计算结果对比------------
                               conditionalPanel("input.Real_PE == 'Yes' ",
                                                tabsetPanel(
                                                  tabPanel("Summary",
                                                           h4(em("不同时间尺度潜在蒸散发计算精度")),
                                                           p("If the slope is less than 1, the model overestimates the PE."),
                                                           tableOutput("Statable")
                                                                
                                                  ),
                                                  tabPanel("Scatter & Linear regression",
                                                           p("If the scatters are under the 1:1 line, the model overestimates the PE."),
                                                           plotOutput("VSplot",height = 400)      
                                                  ),
                                                  
                                                  tabPanel("Process graph", 
                                                           h4(em("")),
                                                           plotOutput("processplot",height = 600)  
                                                  )
                                                )
                                               )  
                        )
                      )      
             )
             #--------待添加功能--------
             
             
  )
)

#-------------------server-----------------
server <- function(input, output) {
  #-----------data input and display--------------
  
  filedata <- reactive({
    req(input$DATA)
    df1 <- read.csv(input$DATA$datapath,
                    header = input$header,
                    sep = input$sep,
                    quote = input$quote)
    if (input$PEmodel == "Penman-Monteith") {
      #年月日；日平均、最高、最低气温（℃）；相对湿度（%）；日照时数（h）；气压（kPa）；风速(m3/s)；站点纬度(°)
      `colnames<-`(df1,c("year","month","day","Ta","Tmax","Tmin","RH","daylight","Pa","Wv"))
    } else if (input$PEmodel == "Hamon"  | input$PEmodel == "Blaney_Criddle" | input$PEmodel == "Kharrufa") {
      #年月日；日平均气温（℃）日照时数（h）
      `colnames<-`(df1,c("year","month","day","Ta","daylight"))
      
    } else if (input$PEmodel == "Linacre" ) {
      #年月日；日平均气温（℃）；相对湿度（%）
      `colnames<-`(df1,c("year","month","day","Ta","RH"))
      
    } else if (input$PEmodel == "Hargreaves") {
      #年月日；日平均、最高、最低气温（℃）；相对湿度（%）；站点纬度(°)
      `colnames<-`(df1,c("year","month","day","Ta","Tmax","Tmin","RH"))
      
    } else if (input$PEmodel == "Romanenko") {
      #年，月；月平均气温（℃）；月平均相对湿度（%）
      `colnames<-`(df1,c("year","month","Ta","RH"))
    }  else if (input$PEmodel == "Thornthwaite") {
      #年,月；月平均气温（℃）；日照时数（h）
      `colnames<-`(df1,c("year","month","Ta","daylight"))
    } else if (input$PEmodel == "Makkink"| input$PEmodel == "FAO-24 Radiation") {
      `colnames<-`(df1,c("year","month","day","Tmax","Tmin","daylight","Pa"))
    } else if ( input$PEmodel == "Priestley-Taylor") {
      `colnames<-`(df1,c("year","month","day","Ta","Tmax","Tmin","daylight","Pa","RH"))
    } else if (input$PEmodel == "Turc") {
      `colnames<-`(df1,c("year","month","day","Tmax","Tmin","daylight","Pa","RH"))
      
    }
    
  })
  output$range <- renderTable({
    rep(input$DATA)
    dat = filedata()
    if (input$Sequencetype == "daily" ) {
      i = paste0(dat[1,1],"-",dat[1,2],"-",dat[1,3])
      m = dim(dat)[1]
      e = paste0(dat[m,1],"-",dat[m,2],"-",dat[m,3])
    } else {
      i = paste0(dat[1,1],"-",dat[1,2])
      m = dim(dat)[1]
      e = paste0(dat[m,1],"-",dat[m,2])
    }

    out = data.frame(i,e)
    colnames(out) = c("From","To")
    out
  })
  output$DATA_show <- DT::renderDataTable({   #upload data
    req(input$DATA)
    df = filedata()
    
    DT::datatable(df)
  })
  
  #---------------simulation-------------------

  datarange <- reactive({
    rep(input$DATA)
    dat = filedata()
    yi = input$period[1]
    ye = input$period[2]
    out = dat[dat[,1] >= yi & dat[,1] <= ye,]
    out
  })

  #----------------PE calculation---------------
  PE_cal <- reactive({
    rep(input$DATA)
    dat = datarange()
    if (input$PEmodel == "Penman-Monteith") {
      lat = input$lat
      J = NULL
      for (j in 1:length(dat[,1])) {
        J[j] = DayofID(dat[j,1],dat[j,2],dat[j,3])
      }
      
      Tmax = dat$Tmax
      Tmin = dat$Tmin
      Ta = (Tmax-Tmin)/2
      RH = dat$RH
      Wv = dat$Wv
      n = dat$daylight
      Pa = dat$Pa
      PE = Penman(Ta,Tmax,Tmin,RH,Wv,Pa,n,J,lat)
      PE[PE<0]=0
    } else if (input$PEmodel == "Hamon" ) {
      D = dat$daylight
      Ta = dat$Ta
      ha = input$ha
      PE = Hamon(D,Ta,ha) 
      
    } else if (input$PEmodel == "Thornthwaite" ) {
      C = input$C
      run = c(31,29,31,30,31,30,31,31,30,31,30,31)
      ping = c(31,28,31,30,31,30,31,31,30,31,30,31)
      
      Ta = dat$Ta
      d = dat$daylight
      PE = NULL
      
      years = as.numeric(names(table(dat[,1])))  #所有的年份
      II = NULL   #计算得到每年的热力指数heat of index：I
      for (t in 1:length(years)) {
        Tm = dat[dat[,1]==years[i],"Ta"]
        II[t] = sum((Tm/5)^1.51)
      }

      for (i in 1:length(dat[,1])) {
        
        y = dat[i,1]
        I = II[years == y]
        if ( (y%%4 == 0) & (y %%100 != 0) |(y %% 400 == 0) ) {   
          days=run
        } else {days=ping}
        N = days[dat[i,2]]   #返回该月的天数
        PE[i] = Thornthwaite(C,I,Ta[i],d[i],N)
      }
      PE[PE < 0] = 0
      
    } else if (input$PEmodel == "Linacre") {
      L = input$L
      lat = input$lat
      Ta = dat$Ta
      RH = dat$RH
      PE = Linacre(Ta,lat,RH,L)
      PE[PE<0] = 0
    } else if (input$PEmodel == "Blaney_Criddle") {
      k = input$k
      Ta = dat$Ta
      yearlight = NULL
      for (i in 1:length(dat[,1])) {
        yearlight[i] = sum(dat[dat[,1] == dat[i,1],"daylight"],na.rm = T)
      }
      p = dat$daylight / yearlight
      PE = Blaney_Criddle(k,Ta,p)
      PE[PE<0] = 0
      
    } else if (input$PEmodel == "Kharrufa") {
      kh = input$kh
      Ta = dat$Ta
      yearlight = NULL
      for (i in 1:length(dat[,1])) {
        yearlight[i] = sum(dat[dat[,1] == dat[i,1],"daylight"])
      }
      p = dat$daylight / yearlight
      PE = Kharrufa(Ta,p,kh)
      PE[PE<0] = 0
      PE[is.na(PE)] = 0
      
    } else if (input$PEmodel == "Hargreaves") {
      k = input$Ha
      lat = input$lat
      Ta = dat$Ta
      RH = dat$RH
      Tmax = dat$Tmax
      Tmin = dat$Tmin
      J = NULL
      for (j in 1:length(dat[,1])) {
        J[j] = DayofID(dat[j,1],dat[j,2],dat[j,3])
      }
      PE = Hargreaves(k,Ta,RH,lat,J,Tmax,Tmin)
      PE[PE<0]=0
    } else if (input$PEmodel == "Romanenko"){
      ro = input$ro
      Ta = dat$Ta
      RH = dat$RH
      PE = Romanenko(Ta,RH,ro)
      PE[PE<0] = 0
    } else if (input$PEmodel == "FAO-24 Radiation") {
      lat = input$lat
      c = input$c
      J = NULL
      for (j in 1:length(dat[,1])) {
        J[j] = DayofID(dat[j,1],dat[j,2],dat[j,3])
      }

      Ta = dat$Ta
      n = dat$daylight
      Pa = dat$Pa
      PE = Penman(Ta,Pa,n,J,lat,c)
      PE[PE<0]=0
      
    } else if (input$PEmodel == "Makkink") {
      lat = input$lat
      
      J = NULL
      for (j in 1:length(dat[,1])) {
        J[j] = DayofID(dat[j,1],dat[j,2],dat[j,3])
      }
      
      Ta = dat$Ta
      n = dat$daylight
      Pa = dat$Pa
      PE = Penman(Ta,Pa,n,J,lat)
      PE[PE<0]=0
      
    } else if (input$PEmodel == "Priestley-Taylor") {
      lat = input$lat
      a = input$alfa
      J = NULL
      for (j in 1:length(dat[,1])) {
        J[j] = DayofID(dat[j,1],dat[j,2],dat[j,3])
      }
      Tmax = dat$Tmax
      Tmin = dat$Tmin
      Ta = (Tmax-Tmin)/2
      RH = dat$RH
      
      n = dat$daylight
      Pa = dat$Pa
      PE = Priestley_Taylor(Ta,Tmax,Tmin,RH,Pa,n,J,lat,a)
      PE[PE<0]=0
      
    } else if (input$PEmodel == "Turc") {
      lat = input$lat
      
      J = NULL
      for (j in 1:length(dat[,1])) {
        J[j] = DayofID(dat[j,1],dat[j,2],dat[j,3])
      }
      Tmax = dat$Tmax
      Tmin = dat$Tmin
      Ta = (Tmax-Tmin)/2
      RH = dat$RH
      
      n = dat$daylight
      Pa = dat$Pa
      PE = Turc(Ta,RH,Pa,n,J,lat)
      PE[PE<0]=0
    }
    
    if (input$PEmodel != "Romanenko" & input$PEmodel != "Thornthwaite") {
      out=data.frame(dat[,c(1:3)],PE)
      `colnames<-`(out,c("y","m","d","PEc"))
    } else {
      out=data.frame(dat[,c(1:2)],PE)
      `colnames<-`(out,c("y","m","PEc"))
    }
     
     
  })
#------reference PE data input---------
  PE_obs <- reactive({
    req(input$RealPE)
    df1 <- read.csv(input$RealPE$datapath,
                    header = input$header,
                    sep = input$sep,
                    quote = input$quote)
    r = input$period
    f = r[1]
    e = r[2]
    out = df1[(df1[,1] >= f & df1[,1] <= e),]
    if (input$Sequencetype2 == "daily") {
      `colnames<-`(out,c("year","month","day","PEr"))
    } else {
      `colnames<-`(out,c("year","month","PEr"))
    }
      
  })
 #-----data rescaled------- 
  monthlydata <- reactive({
    req(input$DATA)
    req(input$RealPE)
    
    PEc = PE_cal()
    PEo = PE_obs()
    
    #****************计算潜在蒸散发整理成月尺度******************
    if (input$PEmodel != "Romanenko" & input$PEmodel != "Thornthwaite") { #daily method
      years = as.numeric(names(table(PEc[,1])))
      pec = NULL
      
      for (i in 1:length(years)) {
        for (j in 1:12) {
          pec[j+(i-1)*12] = sum(PEc[PEc[,1] == years[i] & PEc[,2] == j,4],na.rm = T)
        }
      }
      
    } else {   #当模型为日尺度时；method:Thornthwaite,Romanenko
      pec = PEc[,3]
      YY = PEc[,1]
      MM = PEc[,2]
    }
    #*************观测（observed）蒸散发整理成月尺度***********************
    if (input$Sequencetype2 == "daily") {  #当观测资料输入的为daily尺度时，需要转换
      years = as.numeric(names(table(PEo[,1])))
      peo = NULL
      YY = NULL
      for (i in 1:length(years)) {
        for (j in 1:12) {
          peo[j+(i-1)*12] = sum(PEo[PEo[,1] == years[i] & PEo[,2] == j,4],na.rm = T)
        }
        YY = c(YY,rep(years[i],12))
      }
      
      MM = rep(1:12,length(years))
    } else {
      peo = PEo[,3]
      YY = PEo[,1]
      MM = PEo[,2]
    }
    out = data.frame(YY,MM,pec,peo)
    `colnames<-`(out,c("y","m","pec","peo"))
  })
  
  seasonaldata <- reactive({
    req(input$DATA)
    req(input$RealPE)
    
    dat = monthlydata()
    years = as.numeric(names(table(dat[,1])))
    sc1 = sc2 = sc3 = sc4 = NULL
    so1 = so2 = so3 = so4 = NULL
    
    for (i in 1:length(years)) {
      
      sc1[i] = sum(dat[dat[,1] == years[i] & (dat[,2] >= 3 & dat[,2] <= 5),"pec"],na.rm = T)
      sc2[i] = sum(dat[dat[,1] == years[i] & (dat[,2] >= 6 & dat[,2] <= 8),"pec"],na.rm = T)
      sc3[i] = sum(dat[dat[,1] == years[i] & (dat[,2] >= 9 & dat[,2] <= 11),"pec"],na.rm = T)
      sc4[i] = sum(dat[dat[,1] == years[i] & (dat[,2] == 1 | dat[,2] == 2 | dat[,2] == 12),"pec"],na.rm = T)
      
      so1[i] = sum(dat[dat[,1] == years[i] & (dat[,2] >= 3 & dat[,2] <= 5),"peo"],na.rm = T)
      so2[i] = sum(dat[dat[,1] == years[i] & (dat[,2] >= 6 & dat[,2] <= 8),"peo"],na.rm = T)
      so3[i] = sum(dat[dat[,1] == years[i] & (dat[,2] >= 9 & dat[,2] <= 11),"peo"],na.rm = T)
      so4[i] = sum(dat[dat[,1] == years[i] & (dat[,2] == 1 | dat[,2] == 2 | dat[,2] == 12),"peo"],na.rm = T)
    }
    YY = rep(years,4)
    SS = c(rep("Spring",length(so1)),
           rep("Summer",length(so2)),
           rep("Autumn",length(so3)),
           rep("Winter",length(so4)))
    pec = c(sc1,sc2,sc3,sc4)
    peo = c(so1,so2,so3,so4)
    out = data.frame(YY,SS,pec,peo)
    `colnames<-`(out,c("y","s","pec","peo"))
    
  })
  annualdata <- reactive({
    req(input$DATA)
    req(input$RealPE)

    PEc = PE_cal()
    PEo = PE_obs()
    dat = monthlydata()
    
    years = as.numeric(names(table(dat[,1])))
    pec= NULL
    peo = NULL
    for (i in 1:length(years)) {
      pec[i] = sum(PEc[PEc[,1] == years[i],"PEc"],na.rm = T)
      peo[i] = sum(PEo[PEo[,1] == years[i],"PEr"],na.rm = T)
    }
    out = data.frame(years,pec,peo)
    `colnames<-`(out,c("y","pec","peo"))
    
  })
  output$Statable <- renderTable({
    req(input$DATA)
    req(input$RealPE)
    D = function(P,O) {  #评价指标函数: the index of agreement
      D1 = 1-(sum((P-O)^2,na.rm = T)/sum((abs(P-mean(O,na.rm = T))+abs(O-mean(O,na.rm = T)))^2,na.rm = T)) 
      return(round(D1,digits = 2))
    }
    
    EF = function(P,O) { #evaluation criteria :modeling efficiency
      ef = sum(((O-mean(O,na.rm = T))^2 - (P-O)^2),na.rm = T) / sum((O-mean(O,na.rm = T))^2,na.rm = T)
      return(round(ef,digits = 2))
    }
    
    RMSE = function(P,O) { #均方根误差
      rmse = sqrt(sum((P-O)^2,na.rm = T)/length(P))
      return(round(rmse,digits = 2))
    }
    #------------------------------
    dat_s = seasonaldata()
    lim_s = lm(dat_s[,"peo"]~dat_s[,"pec"])
    intercept_s = lim_s$coefficients[[1]]
    S_s = lim_s$coefficients[[2]]
    Co_s = cor(dat_s[,"peo"],dat_s[,"pec"])
    D_s = D(dat_s[,"pec"],dat_s[,"peo"])
    EF_s = EF(dat_s[,"pec"],dat_s[,"peo"])
    RMSE_s = RMSE(dat_s[,"pec"],dat_s[,"peo"])
    
    dat_m = monthlydata()
    lim_m = lm(dat_m[,"peo"]~dat_m[,"pec"])
    intercept_m = lim_m$coefficients[[1]]
    S_m = lim_m$coefficients[[2]]
    Co_m = cor(dat_m[,"peo"],dat_m[,"pec"])
    D_m = D(dat_m[,"pec"],dat_m[,"peo"])
    EF_m = EF(dat_m[,"pec"],dat_m[,"peo"])
    RMSE_m = RMSE(dat_m[,"pec"],dat_m[,"peo"])
    
    dat_a = annualdata()
    lim_a = lm(dat_a[,"peo"]~dat_a[,"pec"])
    intercept_a = lim_a$coefficients[[1]]
    S_a = lim_a$coefficients[[2]]
    Co_a = cor(dat_a[,"peo"],dat_a[,"pec"])
    D_a = D(dat_a[,"pec"],dat_a[,"peo"])
    EF_a = EF(dat_a[,"pec"],dat_a[,"peo"])
    RMSE_a = RMSE(dat_a[,"pec"],dat_a[,"peo"])
    
    if (input$Sequencetype == "daily" & input$PEmodel != "Romanenko") {
      pe_daily_c = PE_cal()[,4]
      pe_daily_o = PE_obs()[,4]
      lim_d = lm(pe_daily_o ~ pe_daily_c)
      intercept_d = lim_s$coefficients[[1]]
      S_d = lim_d$coefficients[[2]]
      Co_d = cor(pe_daily_o,pe_daily_c)
      D_d = D(pe_daily_c,pe_daily_o)
      EF_d = EF(pe_daily_c,pe_daily_o)
      RMSE_d = RMSE(pe_daily_c,pe_daily_o)
      
      #*******************
      intercept = c(intercept_d,intercept_m,intercept_s,intercept_a)
      S = c(S_d,S_m,S_s,S_a)
      Co = c(Co_d,Co_m,Co_s,Co_a)
      Dd = c(D_d,D_m,D_s,D_a)
      ef = c(EF_d,EF_m,EF_s,EF_a)
      rmse = c(RMSE_d,RMSE_m,RMSE_s,RMSE_a)
      Scale = c("daily","monthly","seasonal","annual")
      out = data.frame(Scale,intercept,S,Co,Dd,ef,rmse)
    } else {
      intercept = c(intercept_m,intercept_s,intercept_a)
      S = c(S_m,S_s,S_a)
      Co = c(Co_m,Co_s,Co_a)
      Dd = c(D_m,D_s,D_a)
      ef = c(EF_m,EF_s,EF_a)
      rmse = c(RMSE_m,RMSE_s,RMSE_a)
      Scale = c("monthly","seasonal","annual")
      out = data.frame(Scale,intercept,S,Co,Dd,ef,rmse)
    }
    `colnames<-`(out,c("Scale","Intercept","Slope","Correlation","index of agreement:D","Modeling efficiency","RMSE"))
  })
  
  output$VSplot <- renderPlot({
    req(input$DATA)
    req(input$RealPE)
    dat_m = monthlydata()
    dat_s = seasonaldata()
    dat_a = annualdata()
    
    opar = par(no.readonly = TRUE)
    par(mfrow = c(1,3))
    #**************monthly******************
    lim_m = c(min(c(dat_m[,"pec"],dat_m[,"peo"])),max(c(dat_m[,"pec"],dat_m[,"peo"])))
    plot(x = dat_m[,"pec"],y = dat_m[,"peo"],main = "monthly",xlim = lim_m,ylim = lim_m,
         xlab = "Calculated",ylab = "Observed")
    ym = xm = min(c(dat_m[,"pec"],dat_m[,"peo"])):max(c(dat_m[,"pec"],dat_m[,"peo"]))
    abline(lm(dat_m[,"peo"]~dat_m[,"pec"]),lwd = 1.5,col = "black")
    points(xm,ym,type = "l",col = "red",lwd = 1.5)
    
    #*************seasonal*********************
    lim_s = c(min(c(dat_s[,"pec"],dat_s[,"peo"])),max(c(dat_s[,"pec"],dat_s[,"peo"])))
    plot(x = dat_s[dat_s[,"s"] == "Spring","pec"],y = dat_s[dat_s[,"s"] == "Spring","peo"],
         main = "seasonal",xlim = lim_s,ylim = lim_s,
         type = "p",pch = 15,col = "gray22",cex=1.5,
         xlab = "Calculated",ylab = "Observed")
    points(x = dat_s[dat_s[,"s"] == "Summer","pec"],y = dat_s[dat_s[,"s"] == "Summer","peo"],
           type = "p",pch = 17,col = "coral4",cex=1.5)
    points(x = dat_s[dat_s[,"s"] == "Autumn","pec"],y = dat_s[dat_s[,"s"] == "Autumn","peo"],
           type = "p",pch = 20,col = "firebrick4",cex=1.5)
    points(x = dat_s[dat_s[,"s"] == "Winter","pec"],y = dat_s[dat_s[,"s"] == "Winter","peo"],
           type = "p",pch = 18,col = "yellow4",cex=1.5)
    legend("topleft",
           c("Spring","Summer","Autumn","Winter"),pch = c(15,17,20,18),cex = c(1.5,1.5,1.5,1.5),
           col = c("gray22","coral4","firebrick4","yellow4"))
    ys = xs = min(c(dat_s[,"pec"],dat_s[,"peo"])):max(c(dat_s[,"pec"],dat_s[,"peo"]))
    
    abline(lm(dat_s[,"peo"]~dat_s[,"pec"]),lwd = 1.5,col = "black")
    points(xs,ys,type = "l",col = "red",lwd = 1.5)
    
    
    #***********annual**********************
    lim_a = c(min(c(dat_a[,"pec"],dat_a[,"peo"])),max(c(dat_a[,"pec"],dat_a[,"peo"])))
    plot(x = dat_a[,"pec"],y = dat_a[,"peo"],
         main = "Annual",xlim = lim_a,ylim = lim_a,
         xlab = "Calculated",ylab = "Observed")
    
    ya = xa = min(c(dat_a[,"pec"],dat_a[,"peo"])):max(c(dat_a[,"pec"],dat_a[,"peo"]))
    abline(lm(dat_a[,"peo"]~dat_a[,"pec"]),lwd = 1.5,col = "black")
    points(xa,ya,type = "l",col = "red",lwd = 1.5)
    
    par(opar)
  })
  output$processplot <- renderPlot({
    req(input$DATA)
    req(input$RealPE)
    dat_m = monthlydata()
    dat_a = annualdata()
    
    opar = par(no.readonly = TRUE)
    par(mfrow = c(2,1))
    x = 1:length(dat_m[,1])
    ylimit = c(0,max(dat_m[,"peo"] ,dat_m[,"pec"] )+10)
    plot(x,dat_m[,"peo"],type = "b",pch = 21,col = "gray11",main = "Monthly",ylim = ylimit,
         xlab = "Calculated",ylab = "Observed")
    points(x,dat_m[,"pec"],type = "l",lty = 2,lwd = 2,col = "red")
    legend("topleft",c("Calculated"),lty = 2,lwd = 2,col = "red")
    
    ylimit2 = c(min(dat_a[,"peo"] ,dat_a[,"pec"] )-100,max(dat_a[,"peo"] ,dat_a[,"pec"] )+100)
    plot(dat_a[,"y"],dat_a[,"peo"],type = "b",pch = 21,col = "gray11",main = "Annual",ylim = ylimit2,
         xlab = "Calculated",ylab = "Observed")
    points(dat_a[,"y"],dat_a[,"pec"],type = "l",lty = 2,lwd = 2,col = "red")
    legend("topleft",c("Calculated"),lty = 2,lwd = 2,col = "red")
    par(opar)
  })



  #----- Downloadable csv ----
  output$downloadData <- downloadHandler(
    filename = function() {
      paste0("潜在蒸散发计算结果.csv")
    },
    content = function(file) {
      write.csv(PE_cal(), file, row.names = FALSE)
    }
  )
}


shinyApp(ui = ui, server = server)
