# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0000156","phosphorelay response regulator activity",0.10002303428214304,0.9858333049307602,-0,"phosphorelay response regulator activity"),
c("GO:0004970","ionotropic glutamate receptor activity",0.04584981516793261,0.932184002941164,0.52275345,"phosphorelay response regulator activity"),
c("GO:0016989","sigma factor antagonist activity",0.008282388486621413,0.9752982456395392,0.15541818,"phosphorelay response regulator activity"),
c("GO:0070038","rRNA (pseudouridine-N3-)-methyltransferase activity",0.016701495172097586,0.8847892426901462,0.19147382,"phosphorelay response regulator activity"),
c("GO:0004556","alpha-amylase activity",0.02656161167350294,0.9334924925703809,0.0167541,"alpha-amylase activity"),
c("GO:0003964","RNA-directed DNA polymerase activity",0.1316158756735012,0.8445825203446364,0.47804567,"alpha-amylase activity"),
c("GO:0004035","alkaline phosphatase activity",0.019337422046017377,0.9222857548825663,0.29700496,"alpha-amylase activity"),
c("GO:0004061","arylformamidase activity",0.03657758692160272,0.9297738469409424,0.61099483,"alpha-amylase activity"),
c("GO:0004519","endonuclease activity",1.2573644624995106,0.9011705022369081,0.43995601,"alpha-amylase activity"),
c("GO:0008233","peptidase activity",4.315214635541,0.903831374432116,0.21714218,"alpha-amylase activity"),
c("GO:0008253","5'-nucleotidase activity",0.05164119807142027,0.9179802152199685,0.53403136,"alpha-amylase activity"),
c("GO:0008725","DNA-3-methyladenine glycosylase activity",0.02719051538823484,0.9120233893495527,0.46607468,"alpha-amylase activity"),
c("GO:0008859","exoribonuclease II activity",0.01879328361457543,0.894838877698871,0.58497104,"alpha-amylase activity"),
c("GO:0009002","serine-type D-Ala-D-Ala carboxypeptidase activity",0.09146174066985793,0.9236359097770649,0.38744842,"alpha-amylase activity"),
c("GO:0009039","urease activity",0.032276432385631944,0.9302331021593548,0.14921291,"alpha-amylase activity"),
c("GO:0016920","pyroglutamyl-peptidase activity",0.010524566947839492,0.932480696082662,0.31935151,"alpha-amylase activity"),
c("GO:0047746","chlorophyllase activity",0.0024964743110879335,0.9332605715362645,0.12893285,"alpha-amylase activity"),
c("GO:0015625","ABC-type ferric hydroxamate transporter activity",5.468727954190434E-06,0.9491163838503254,0.00342306,"ABC-type ferric hydroxamate transporter activity"),
c("GO:0005381","iron ion transmembrane transporter activity",0.05885991897095165,0.9326136981831498,0.48319221,"ABC-type ferric hydroxamate transporter activity"),
c("GO:0015116","sulfate transmembrane transporter activity",0.04804550944154006,0.9264327924808161,0.29922617,"ABC-type ferric hydroxamate transporter activity"),
c("GO:0015137","citrate transmembrane transporter activity",0.004128889605413778,0.9329050968071377,0.5643334,"ABC-type ferric hydroxamate transporter activity"),
c("GO:0015386","potassium:proton antiporter activity",0.007360907826340325,0.9217446846561459,0.52211226,"ABC-type ferric hydroxamate transporter activity"),
c("GO:0015416","ABC-type phosphonate transporter activity",0.0218393650850595,0.9381544695133303,0.3653369,"ABC-type ferric hydroxamate transporter activity"),
c("GO:0015558","secondary active p-aminobenzoyl-glutamate transmembrane transporter activity",0.005542555781572005,0.9329851700962375,0.5215295,"ABC-type ferric hydroxamate transporter activity"),
c("GO:0043858","arginine:ornithine antiporter activity",0.0008066373732430891,0.9366575488953609,0.24187631,"ABC-type ferric hydroxamate transporter activity"),
c("GO:0016151","nickel cation binding",0.11700070021592726,0.9952940027738303,0.03657963,"nickel cation binding"),
c("GO:0030975","thiamine binding",0.018580003224362002,0.9946895023993112,0.1802912,"nickel cation binding"),
c("GO:0019808","polyamine binding",0.023865528792087056,0.9982111508052154,0,"polyamine binding"),
c("GO:0030983","mismatched DNA binding",0.11158392517730163,0.9973303539641457,0.03273832,"mismatched DNA binding"),
c("GO:0043800","hexulose-6-phosphate isomerase activity",0.0004019515046329969,0.9620751724280723,0.01391246,"hexulose-6-phosphate isomerase activity"),
c("GO:0008736","L-fucose isomerase activity",0.00534021284726696,0.9538400802832737,0.54782691,"hexulose-6-phosphate isomerase activity"),
c("GO:0050044","galactose-6-phosphate isomerase activity",0.0007601531856324704,0.9571687987388039,0.3001879,"hexulose-6-phosphate isomerase activity"),
c("GO:0047451","3-hydroxyoctanoyl-[acyl-carrier-protein] dehydratase activity",8.749964726704695E-05,0.9139542193302933,0.01310388,"3-hydroxyoctanoyl-[acyl-carrier-protein] dehydratase activity"),
c("GO:0004049","anthranilate synthase activity",0.02841004172201931,0.8827820222351414,0.49453414,"3-hydroxyoctanoyl-[acyl-carrier-protein] dehydratase activity"),
c("GO:0004121","cystathionine beta-lyase activity",0.011621046902654674,0.9006106600800325,0.28126437,"3-hydroxyoctanoyl-[acyl-carrier-protein] dehydratase activity"),
c("GO:0004160","dihydroxy-acid dehydratase activity",0.03337838106840132,0.8867680667243109,0.57087716,"3-hydroxyoctanoyl-[acyl-carrier-protein] dehydratase activity"),
c("GO:0004612","phosphoenolpyruvate carboxykinase (ATP) activity",0.01797570878542396,0.8826932848647161,0.54342705,"3-hydroxyoctanoyl-[acyl-carrier-protein] dehydratase activity"),
c("GO:0008659","(3R)-hydroxymyristoyl-[acyl-carrier-protein] dehydratase activity",0.028289729707027117,0.8877347756146521,0.53922257,"3-hydroxyoctanoyl-[acyl-carrier-protein] dehydratase activity"),
c("GO:0016832","aldehyde-lyase activity",0.18000591497615526,0.8705750725028275,0.62985729,"3-hydroxyoctanoyl-[acyl-carrier-protein] dehydratase activity"),
c("GO:0016840","carbon-nitrogen lyase activity",0.26838876180780397,0.8819034334954188,0.43506717,"3-hydroxyoctanoyl-[acyl-carrier-protein] dehydratase activity"),
c("GO:0047434","indolepyruvate decarboxylase activity",0.0005304666115564722,0.900526022174756,0.38952912,"3-hydroxyoctanoyl-[acyl-carrier-protein] dehydratase activity"),
c("GO:0047456","2-methylisocitrate dehydratase activity",0.012914401063820712,0.8921255544302045,0.40317592,"3-hydroxyoctanoyl-[acyl-carrier-protein] dehydratase activity"),
c("GO:0047605","acetolactate decarboxylase activity",0.0034179549713690216,0.8917817652769138,0.45549652,"3-hydroxyoctanoyl-[acyl-carrier-protein] dehydratase activity"),
c("GO:0050150","o-pyrocatechuate decarboxylase activity",7.656219135866609E-05,0.9082772266087753,0.2308698,"3-hydroxyoctanoyl-[acyl-carrier-protein] dehydratase activity"),
c("GO:0047482","UDP-N-acetylmuramoyl-L-alanyl-D-glutamate-L-lysine ligase activity",0.00046210751212909175,0.9682611334467286,-0,"UDP-N-acetylmuramoyl-L-alanyl-D-glutamate-L-lysine ligase activity"),
c("GO:0047760","butyrate-CoA ligase activity",0.0005085916997397104,0.9681382921113971,0.26782947,"UDP-N-acetylmuramoyl-L-alanyl-D-glutamate-L-lysine ligase activity"),
c("GO:0047810","D-alanine:2-oxoglutarate aminotransferase activity",0.002518349222904695,0.909761162675653,0.01826331,"D-alanine:2-oxoglutarate aminotransferase activity"),
c("GO:0000310","xanthine phosphoribosyltransferase activity",0.009075354039979026,0.8899571019438485,0.14814038,"D-alanine:2-oxoglutarate aminotransferase activity"),
c("GO:0003961","O-acetylhomoserine aminocarboxypropyltransferase activity",0.006390208614471524,0.8954998217632519,0.14504001,"D-alanine:2-oxoglutarate aminotransferase activity"),
c("GO:0004020","adenylylsulfate kinase activity",0.026591689677250986,0.8816857993029457,0.32275728,"D-alanine:2-oxoglutarate aminotransferase activity"),
c("GO:0004314","[acyl-carrier-protein] S-malonyltransferase activity",0.029036211072774112,0.8571750750184456,0.44212262,"D-alanine:2-oxoglutarate aminotransferase activity"),
c("GO:0004810","CCA tRNA nucleotidyltransferase activity",0.019643670811452043,0.8593167592125056,0.68968531,"D-alanine:2-oxoglutarate aminotransferase activity"),
c("GO:0004849","uridine kinase activity",0.01855812831254524,0.8777240651372042,0.28438204,"D-alanine:2-oxoglutarate aminotransferase activity"),
c("GO:0004850","uridine phosphorylase activity",0.011295657589380342,0.8887444906387912,0.56366403,"D-alanine:2-oxoglutarate aminotransferase activity"),
c("GO:0008410","CoA-transferase activity",0.08223326224716157,0.8988361729678626,0.13740658,"D-alanine:2-oxoglutarate aminotransferase activity"),
c("GO:0008728","GTP diphosphokinase activity",0.012053076411035719,0.8897052948070645,0.24263968,"D-alanine:2-oxoglutarate aminotransferase activity"),
c("GO:0008808","cardiolipin synthase activity",0.016868291374700396,0.8876824333648726,0.25620408,"D-alanine:2-oxoglutarate aminotransferase activity"),
c("GO:0008865","fructokinase activity",0.005091385725351295,0.8899230236744382,0.29548467,"D-alanine:2-oxoglutarate aminotransferase activity"),
c("GO:0008955","peptidoglycan glycosyltransferase activity",0.0872453514171771,0.8824496508788409,0.44529653,"D-alanine:2-oxoglutarate aminotransferase activity"),
c("GO:0008959","phosphate acetyltransferase activity",0.01001050652014559,0.8730635512360723,0.55669156,"D-alanine:2-oxoglutarate aminotransferase activity"),
c("GO:0009024","tagatose-6-phosphate kinase activity",0.006267162235502239,0.8904347754851468,0.14487196,"D-alanine:2-oxoglutarate aminotransferase activity"),
c("GO:0009032","thymidine phosphorylase activity",0.012610886662363142,0.8881251002232389,0.57392837,"D-alanine:2-oxoglutarate aminotransferase activity"),
c("GO:0016775","phosphotransferase activity, nitrogenous group as acceptor",0.9908952242036275,0.856263783066161,0.35321717,"D-alanine:2-oxoglutarate aminotransferase activity"),
c("GO:0034386","4-aminobutyrate:2-oxoglutarate transaminase activity",0.008465590873086794,0.9043675049536468,0.49607265,"D-alanine:2-oxoglutarate aminotransferase activity"),
c("GO:0043772","acyl-phosphate glycerol-3-phosphate acyltransferase activity",0.01800578678917201,0.8758563322530445,0.39105518,"D-alanine:2-oxoglutarate aminotransferase activity"),
c("GO:0047200","tetrahydrodipicolinate N-acetyltransferase activity",0.003251158768766213,0.8799003061548657,0.1394192,"D-alanine:2-oxoglutarate aminotransferase activity"),
c("GO:0047348","glycerol-3-phosphate cytidylyltransferase activity",0.001998820067256604,0.8878168717314294,0.22137094,"D-alanine:2-oxoglutarate aminotransferase activity"),
c("GO:0051996","squalene synthase activity",0.006321849515044143,0.8942918451570331,0.48780104,"D-alanine:2-oxoglutarate aminotransferase activity"),
c("GO:0097163","sulfur carrier activity",0.015285094631962266,0.9862973693287981,-0,"sulfur carrier activity"),
c("GO:0102223","4,4'-diapophytoene desaturase (4,4'-diaponeurosporene-forming)",0.0003499985890681878,0.8792351038511868,0.01383497,"4,4'-diapophytoene desaturase (4,4'-diaponeurosporene-forming)"),
c("GO:0003857","3-hydroxyacyl-CoA dehydrogenase activity",0.030731516738573147,0.8120526800766488,0.50265024,"4,4'-diapophytoene desaturase (4,4'-diaponeurosporene-forming)"),
c("GO:0003858","3-hydroxybutyrate dehydrogenase activity",0.010694097514419396,0.8221446927808586,0.46168537,"4,4'-diapophytoene desaturase (4,4'-diaponeurosporene-forming)"),
c("GO:0003863","3-methyl-2-oxobutanoate dehydrogenase (2-methylpropanoyl-transferring) activity",0.010770659705778061,0.8234159182868491,0.46720441,"4,4'-diapophytoene desaturase (4,4'-diaponeurosporene-forming)"),
c("GO:0003960","NADPH:quinone reductase activity",0.006778488299219044,0.8543032862701752,0.19798063,"4,4'-diapophytoene desaturase (4,4'-diaponeurosporene-forming)"),
c("GO:0004029","aldehyde dehydrogenase (NAD+) activity",0.03624126015242001,0.8046394374158732,0.64056173,"4,4'-diapophytoene desaturase (4,4'-diaponeurosporene-forming)"),
c("GO:0004322","ferroxidase activity",0.01969562372701685,0.8526373419753304,0.21938594,"4,4'-diapophytoene desaturase (4,4'-diaponeurosporene-forming)"),
c("GO:0004324","ferredoxin-NADP+ reductase activity",0.02174366234586117,0.8518390355010027,0.22311892,"4,4'-diapophytoene desaturase (4,4'-diaponeurosporene-forming)"),
c("GO:0004355","glutamate synthase (NADPH) activity",0.016291340575533306,0.8541453777384981,0.21581266,"4,4'-diapophytoene desaturase (4,4'-diaponeurosporene-forming)"),
c("GO:0004365","glyceraldehyde-3-phosphate dehydrogenase (NAD+) (phosphorylating) activity",0.019897966661321896,0.8154036665025155,0.5958776,"4,4'-diapophytoene desaturase (4,4'-diaponeurosporene-forming)"),
c("GO:0004392","heme oxygenase (decyclizing) activity",0.014667128373138747,0.8503242264770704,0.47279729,"4,4'-diapophytoene desaturase (4,4'-diaponeurosporene-forming)"),
c("GO:0004420","hydroxymethylglutaryl-CoA reductase (NADPH) activity",0.01595774817032769,0.818441489691389,0.47469971,"4,4'-diapophytoene desaturase (4,4'-diaponeurosporene-forming)"),
c("GO:0004497","monooxygenase activity",1.2376879793203333,0.8103345239595742,0.3352805,"4,4'-diapophytoene desaturase (4,4'-diaponeurosporene-forming)"),
c("GO:0004517","nitric-oxide synthase activity",0.007328095458615183,0.8555189184123552,0.19090178,"4,4'-diapophytoene desaturase (4,4'-diaponeurosporene-forming)"),
c("GO:0004604","phosphoadenylyl-sulfate reductase (thioredoxin) activity",0.018402269565850813,0.8465590920357022,0.59143066,"4,4'-diapophytoene desaturase (4,4'-diaponeurosporene-forming)"),
c("GO:0004783","sulfite reductase (NADPH) activity",0.008599574707964458,0.852363616264947,0.20822786,"4,4'-diapophytoene desaturase (4,4'-diaponeurosporene-forming)"),
c("GO:0008106","alcohol dehydrogenase (NADP+) activity",0.010300349101717683,0.8224844081943441,0.43562948,"4,4'-diapophytoene desaturase (4,4'-diaponeurosporene-forming)"),
c("GO:0008752","FMN reductase activity",0.003335924052056165,0.8580995456049708,0.16551428,"4,4'-diapophytoene desaturase (4,4'-diaponeurosporene-forming)"),
c("GO:0008812","choline dehydrogenase activity",0.007812077882561036,0.826888189720319,0.20713358,"4,4'-diapophytoene desaturase (4,4'-diaponeurosporene-forming)"),
c("GO:0008863","formate dehydrogenase (NAD+) activity",0.028508478825194736,0.8119005599257934,0.61754299,"4,4'-diapophytoene desaturase (4,4'-diaponeurosporene-forming)"),
c("GO:0008924","malate dehydrogenase (quinone) activity",0.006961690685684424,0.8228495515186107,0.42666141,"4,4'-diapophytoene desaturase (4,4'-diaponeurosporene-forming)"),
c("GO:0016652","oxidoreductase activity, acting on NAD(P)H, NAD(P) as acceptor",0.04205725233170154,0.8399113927606824,0.58619276,"4,4'-diapophytoene desaturase (4,4'-diaponeurosporene-forming)"),
c("GO:0016661","oxidoreductase activity, acting on other nitrogenous compounds as donors",0.10092810875856158,0.8382802850122923,0.24689769,"4,4'-diapophytoene desaturase (4,4'-diaponeurosporene-forming)"),
c("GO:0043878","glyceraldehyde-3-phosphate dehydrogenase (NAD+) (non-phosphorylating) activity",0.015542124845809214,0.8177410735739588,0.54243941,"4,4'-diapophytoene desaturase (4,4'-diaponeurosporene-forming)"),
c("GO:0046857","oxidoreductase activity, acting on other nitrogenous compounds as donors, with NAD or NADP as acceptor",0.04765176102883835,0.8249233686283589,0.26021372,"4,4'-diapophytoene desaturase (4,4'-diaponeurosporene-forming)"),
c("GO:0047112","pyruvate oxidase activity",0.0008695277447162791,0.8434624917888623,0.1723729,"4,4'-diapophytoene desaturase (4,4'-diaponeurosporene-forming)"),
c("GO:0047936","glucose 1-dehydrogenase [NAD(P)] activity",0.011407766512441247,0.8215567566007147,0.46427269,"4,4'-diapophytoene desaturase (4,4'-diaponeurosporene-forming)"),
c("GO:0052693","epoxyqueuosine reductase activity",0.014801112208016413,0.8548963044343075,0.20641124,"4,4'-diapophytoene desaturase (4,4'-diaponeurosporene-forming)"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="revigo_treemap.pdf", width=16, height=9 ) # width and height are in inches

# check the tmPlot command documentation for all possible parameters - there are a lot more
treemap(
  stuff,
  index = c("representative","description"),
  vSize = "uniqueness",
  type = "categorical",
  vColor = "representative",
  title = "Revigo TreeMap",
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
								 # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)

dev.off()

