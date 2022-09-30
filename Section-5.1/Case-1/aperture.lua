local aperture= { 
-- idEl, gpt, openValue 
{   15019, 1, 1.264966326E-03}, 
{   15019, 2, 1.257073716E-03}, 
{   15020, 1, 1.257031108E-03}, 
{   15020, 2, 1.275764196E-03}, 
{   15021, 1, 1.275714952E-03}, 
{   15021, 2, 1.265268307E-03}, 
{   15022, 1, 1.265234780E-03}, 
{   15022, 2, 1.287579653E-03}, 
{   15631, 1, 2.501088020E-04}, 
{   15631, 2, 1.521000831E-05}, 
{   15632, 1, 1.556111238E-05}, 
{   15632, 2, 9.992934793E-06}, 
{   16565, 1, 1.002428689E-05}, 
{   16565, 2, 9.999856047E-06}, 
{   16566, 1, 8.207364008E-04}, 
{   16566, 2, 8.801803924E-04}, 
{   16568, 1, 8.762950310E-04}, 
{   16568, 2, 8.207474020E-04}, 
{   17951, 1, 8.167758351E-04}, 
{   17951, 2, 7.877525641E-04}, 
{   17954, 1, 9.416986140E-04}, 
{   17954, 2, 8.167803753E-04}, 
{   18475, 1, 1.001035525E-05}, 
{   18475, 2, 1.003472789E-05}, 
{   19606, 1, 1.106343349E-03}, 
{   19606, 2, 1.156909624E-03}, 
{   19610, 1, 1.103641582E-03}, 
{   19610, 2, 1.106363605E-03}, 
{   20247, 1, 1.000783504E-05}, 
{   20247, 2, 1.002837234E-05}, 
{   20896, 1, 9.669294232E-04}, 
{   20896, 2, 1.047644299E-03}, 
{   20899, 1, 9.651139262E-04}, 
{   20899, 2, 9.669990395E-04}, 
{   21495, 1, 8.002316463E-04}, 
{   21495, 2, 8.229900268E-04}, 
{   21498, 1, 8.581678267E-04}, 
{   21498, 2, 8.002070244E-04}, 
{   21499, 1, 1.000529119E-05}, 
{   21499, 2, 1.002362023E-05}, 
{   21636, 1, 1.003308171E-05}, 
{   21636, 2, 1.000066277E-05}, 
{   21692, 1, 1.003053239E-05}, 
{   21692, 2, 1.000812244E-05}, 
{   21693, 1, 1.102183363E-03}, 
{   21693, 2, 1.057488727E-03}, 
{   21698, 1, 1.057551708E-03}, 
{   21698, 2, 1.086640172E-03}, 
{   22170, 1, 1.183327171E-03}, 
{   22170, 2, 1.197524020E-03}, 
{   22173, 1, 1.207740046E-03}, 
{   22173, 2, 1.183332293E-03}, 
{   22979, 1, 1.000698376E-05}, 
{   22979, 2, 1.002938370E-05}, 
{   23664, 1, 1.149448683E-03}, 
{   23664, 2, 1.195133897E-03}, 
{   23667, 1, 1.156605431E-03}, 
{   23667, 2, 1.149457297E-03}, 
{   23693, 1, 1.025472375E-05}, 
{   23693, 2, 9.998632777E-06}, 
{   23695, 1, 8.276684821E-05}, 
{   23695, 2, 1.010133019E-05}, 
{   24410, 1, 1.003081161E-05}, 
{   24410, 2, 1.001607507E-05}, 
{   27831, 1, 1.002258432E-05}, 
{   27831, 2, 1.000392695E-05}, 
{   27881, 1, 1.000432712E-05}, 
{   27881, 2, 1.002960107E-05}, 
{   29097, 1, 1.000749307E-05}, 
{   29097, 2, 1.003272700E-05}, 
{   29214, 1, 4.666039022E-04}, 
{   29214, 2, 3.677973582E-04}, 
{   29216, 1, 6.182196666E-04}, 
{   29216, 2, 4.665719171E-04}, 
{   29260, 1, 1.000476095E-05}, 
{   29260, 2, 4.847028322E-05}, 
{   29265, 1, 4.826229269E-05}, 
{   29265, 2, 3.880609293E-04}, 
{   29302, 1, 6.633053417E-04}, 
{   29302, 2, 6.250616279E-04}, 
{   29304, 1, 7.846626686E-04}, 
{   29304, 2, 6.632995792E-04}, 
{   29488, 1, 1.033859211E-03}, 
{   29488, 2, 1.113552367E-03}, 
{   29491, 1, 1.059748582E-03}, 
{   29491, 2, 1.033845940E-03}, 
{   31358, 1, 1.137860469E-03}, 
{   31358, 2, 1.157637569E-03}, 
{   31362, 1, 1.147440751E-03}, 
{   31362, 2, 1.137878979E-03}, 
{   33944, 1, 1.003223770E-05}, 
{   33944, 2, 1.000473731E-05}, 
{   35818, 1, 1.178948325E-03}, 
{   35818, 2, 1.157113118E-03}, 
{   35822, 1, 1.002624595E-05}, 
{   35822, 2, 1.001185410E-05}, 
{   35823, 1, 1.157102990E-03}, 
{   35823, 2, 1.180068590E-03}, 
{   36396, 1, 1.143875765E-03}, 
{   36396, 2, 1.165633206E-03}, 
{   36400, 1, 1.161921420E-03}, 
{   36400, 2, 1.143911388E-03}, 
{   37090, 1, 1.003428042E-05}, 
{   37090, 2, 1.000580232E-05}, 
{   37391, 1, 1.049736864E-03}, 
{   37391, 2, 1.104459167E-03}, 
{   37394, 1, 1.043060794E-03}, 
{   37394, 2, 1.049751882E-03}, 
{   37621, 1, 7.371864631E-04}, 
{   37621, 2, 6.837813999E-04}, 
{   37623, 1, 8.392485906E-04}, 
{   37623, 2, 7.370105013E-04}, 
{   38093, 1, 1.000407428E-05}, 
{   38093, 2, 1.001808232E-05}, 
{   38200, 1, 1.135543105E-03}, 
{   38200, 2, 1.100485912E-03}, 
{   38202, 1, 1.100485213E-03}, 
{   38202, 2, 1.117900247E-03}, 
{   40254, 1, 1.154526835E-03}, 
{   40254, 2, 1.179168350E-03}, 
{   40255, 1, 1.000505108E-05}, 
{   40255, 2, 1.002415229E-05}, 
{   40256, 1, 1.174785895E-03}, 
{   40256, 2, 1.154551865E-03}, 
{   40360, 1, 1.002817135E-05}, 
{   40360, 2, 1.000610246E-05}, 
{   40630, 1, 1.001757209E-05}, 
{   40630, 2, 1.003298530E-05}, 
{   41332, 1, 1.000487737E-05}, 
{   41332, 2, 1.003230045E-05}, 
{   41773, 1, 9.855406824E-04}, 
{   41773, 2, 1.042305958E-03}, 
{   41776, 1, 9.804121219E-04}, 
{   41776, 2, 9.855222888E-04}, 
{   41795, 1, 5.720140762E-04}, 
{   41795, 2, 4.966406850E-04}, 
{   41797, 1, 6.941099418E-04}, 
{   41797, 2, 5.718041793E-04}, 
{   41798, 1, 1.000896464E-05}, 
{   41798, 2, 1.003529269E-05}, 
{   42324, 1, 1.003168109E-05}, 
{   42324, 2, 1.000689736E-05}, 
{   43235, 1, 8.767505642E-04}, 
{   43235, 2, 9.618708864E-04}, 
{   43238, 1, 8.485685103E-04}, 
{   43238, 2, 8.768847911E-04}, 
{   43247, 1, 3.634972381E-04}, 
{   43247, 2, 4.905918031E-04}, 
{   43248, 1, 9.998117093E-06}, 
{   43248, 2, 1.003081616E-05}, 
{   43249, 1, 2.567116753E-04}, 
{   43249, 2, 3.638389462E-04}, 
{   44189, 1, 1.019742831E-05}, 
{   44189, 2, 1.000589691E-05}, 
} 
return aperture