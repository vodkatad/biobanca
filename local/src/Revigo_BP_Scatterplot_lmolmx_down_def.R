# A plotting R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# --------------------------------------------------------------------------
# If you don't have the ggplot2 package installed, uncomment the following line:
# install.packages( "ggplot2" );
library( ggplot2 );

# --------------------------------------------------------------------------
# If you don't have the scales package installed, uncomment the following line:
# install.packages( "scales" );
library( scales );

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","plot_X","plot_Y","log_size","value","uniqueness","dispensability");
revigo.data <- rbind(c("GO:0000184","nuclear-transcribed mRNA catabolic process, nonsense-mediated decay",0.21958223073025168,4.89774938749947,-4.077866106991881,1.6020599913279623,-2.039002545851808,0.9231081965188392,0.30239641),
c("GO:0002385","mucosal immune response",0.18017003547097574,4.982762566938707,5.4897182908867075,1.5185139398778875,-1.3205059544447137,0.8787882210962535,0.35962382),
c("GO:0002576","platelet degranulation",0.045042508867743934,-5.187017594006661,5.501722688093441,0.9542425094393249,-1.5038411989165439,0.9396011717163293,0.30340512),
c("GO:0003018","vascular process in circulatory system",1.4751421654186139,-0.5545126947573357,1.8595230811742642,2.419955748489758,-1.332900306718457,0.9931775458834168,-0),
c("GO:0006614","SRP-dependent cotranslational protein targeting to membrane",0.09571533134395585,-4.318222795478114,5.490757000240016,1.255272505103306,-2.5400081268795915,0.8781466850958444,0.36724605),
c("GO:0006882","intracellular zinc ion homeostasis",0.1351275266032318,1.535452775547182,8.449723365570573,1.3979400086720377,-1.5167532328761149,1,-0),
c("GO:0006959","humoral immune response",1.2274083666460223,4.49937335783016,5.0919849572285605,2.3404441148401185,-1.9919139461923692,0.8583785596356128,0.42717639),
c("GO:0007187","G protein-coupled receptor signaling pathway, coupled to cyclic nucleotide second messenger",0.31529756207420756,-1.2054472354554495,8.48812573438985,1.7558748556724915,-1.3205059544447137,0.9309701816555498,0.0729565),
c("GO:0010273","detoxification of copper ion",0.016890940825403974,6.640187354088043,0.019678235765588217,0.6020599913279624,-2.455568713867886,0.8546173812089718,0.18570112),
c("GO:0010466","negative regulation of peptidase activity",1.2217780530375542,5.412939656622149,-3.643640330429071,2.3384564936046046,-2.910829389601365,0.8516390663217286,0.02781909),
c("GO:0010742","macrophage derived foam cell differentiation",0.045042508867743934,1.5888239830641808,-1.416144946459828,0.9542425094393249,-1.8018375900809653,0.9887073517460372,0.00329212),
c("GO:0010743","regulation of macrophage derived foam cell differentiation",0.17453972186250774,1.3034224409941737,-5.934476913290593,1.505149978319906,-1.962765637860222,0.9739162509081791,0.02327061),
c("GO:0010744","positive regulation of macrophage derived foam cell differentiation",0.10134564495242385,1.9075304277080551,-5.779989705270174,1.2787536009528289,-1.8778424273270073,0.9746260111662896,0.33732107),
c("GO:0030193","regulation of blood coagulation",0.38849163898429145,-3.13632681183774,-5.06037515316243,1.845098040014257,-1.9963533565816616,0.7546954953429221,-0),
c("GO:0030198","extracellular matrix organization",1.5427059287202296,-5.901547050970343,0.9687811848002582,2.439332693830263,-2.189488090816805,0.9247819699090056,0.19398109),
c("GO:0034340","response to type I interferon",0.29840662124880357,5.466413939847427,3.4899926655426357,1.7323937598229686,-2.039002545851808,0.7788103108265098,0),
c("GO:0034367","protein-containing complex remodeling",0.17453972186250774,-6.2761266432381015,-1.1807779452006595,1.505149978319906,-1.332900306718457,0.940883712960367,0.3059692),
c("GO:0034369","plasma lipoprotein particle remodeling",0.16327909464557178,-5.334203044978732,-0.8581650100025032,1.4771212547196624,-1.440980684737327,0.9025223355334084,0.13096816),
c("GO:0035455","response to interferon-alpha",0.10697595856089183,6.595738149502522,1.9220337825251497,1.3010299956639813,-1.3205059544447137,0.8813375545233568,0.46531119),
c("GO:0042073","intraciliary transport",0.26462473959799565,-5.244950880937949,3.8942175003995705,1.6812412373755872,-3.081307724560586,0.8053286772367826,-0),
c("GO:0043062","extracellular structure organization",1.5483362423286977,-7.022049693485292,0.8584679856646509,2.4409090820652177,-2.75439999142571,0.9417075896664499,0.16012992),
c("GO:0050818","regulation of coagulation",0.41101289341816344,-4.425625014990583,-4.531152191310973,1.8692317197309762,-1.837334748937521,0.9542655285431618,0.2710095),
c("GO:0070972","protein localization to endoplasmic reticulum",0.37723101176735546,-4.149552235785949,5.993951809675731,1.8325089127062364,-1.5522959271930301,0.9178539441886506,0.49157458),
c("GO:0090066","regulation of anatomical structure size",2.792635549800124,-2.0714090808235204,-5.578413141253178,2.696356388733332,-1.332900306718457,0.9377952297603857,0.33290402),
c("GO:0090077","foam cell differentiation",0.05067282247621192,-0.03209796567553343,-1.8509650330179652,1,-1.8018375900809653,0.9886559495241152,0.16909815),
c("GO:0097501","stress response to metal ion",0.05067282247621192,6.235054560393486,0.2861675089325503,1,-2.039002545851808,0.8718024528156346,0.46621021),
c("GO:0097529","myeloid leukocyte migration",0.7826135915770508,4.247624045201558,5.695017716907694,2.146128035678238,-2.0544389160796306,0.8840543688752307,0.33941762),
c("GO:0098742","cell-cell adhesion via plasma-membrane adhesion molecules",1.4695118518101458,0.09106686752214799,6.068345702444863,2.4183012913197452,-1.9919139461923692,0.9726623584734088,0.00442195));

one.data <- data.frame(revigo.data);
names(one.data) <- revigo.names;
one.data <- one.data [(one.data$plot_X != "null" & one.data$plot_Y != "null"), ];
one.data$plot_X <- as.numeric( as.character(one.data$plot_X) );
one.data$plot_Y <- as.numeric( as.character(one.data$plot_Y) );
one.data$log_size <- as.numeric( as.character(one.data$log_size) );
one.data$value <- as.numeric( as.character(one.data$value) );
one.data$frequency <- as.numeric( as.character(one.data$frequency) );
one.data$uniqueness <- as.numeric( as.character(one.data$uniqueness) );
one.data$dispensability <- as.numeric( as.character(one.data$dispensability) );
#head(one.data);


# --------------------------------------------------------------------------
# Names of the axes, sizes of the numbers and letters, names of the columns,
# etc. can be changed below

p1 <- ggplot( data = one.data );
p1 <- p1 + geom_point( aes( plot_X, plot_Y, colour = value, size = log_size), alpha = I(0.6) );
p1 <- p1 + scale_colour_gradientn( colours = c("blue", "green", "yellow", "red"), limits = c( min(one.data$value), 0) );
p1 <- p1 + geom_point( aes(plot_X, plot_Y, size = log_size), shape = 21, fill = "transparent", colour = I (alpha ("black", 0.6) ));
p1 <- p1 + scale_size( range=c(5, 30)) + theme_bw(); # + scale_fill_gradientn(colours = heat_hcl(7), limits = c(-300, 0) );
ex <- one.data [ one.data$dispensability < 0.15, ];
p1 <- p1 + geom_text( data = ex, aes(plot_X, plot_Y, label = description), colour = I(alpha("black", 0.85)), size = 3 );
p1 <- p1 + labs (y = "semantic space x", x = "semantic space y");
p1 <- p1 + theme(legend.key = element_blank()) ;
one.x_range = max(one.data$plot_X) - min(one.data$plot_X);
one.y_range = max(one.data$plot_Y) - min(one.data$plot_Y);
p1 <- p1 + xlim(min(one.data$plot_X)-one.x_range/10,max(one.data$plot_X)+one.x_range/10);
p1 <- p1 + ylim(min(one.data$plot_Y)-one.y_range/10,max(one.data$plot_Y)+one.y_range/10);


# --------------------------------------------------------------------------
# Output the plot to screen

p1;

# Uncomment the line below to also save the plot to a file.
# The file type depends on the extension (default=pdf).

# ggsave("/path_to_your_file/revigo-plot.pdf");

