#include "phys_const.hh"

GPHYS_CONST::GPHYS_CONST()
{;}

GPHYS_CONST::~GPHYS_CONST()
{;}

tuple<ld1, array<ld1, 2*nmax+1>> GPHYS_CONST::cpotss(G4int sigen, G4int k0, G4int l0, G4int n0, G4String crystaltype)
{
    array<ld1, 2*nmax+1> cpot;
    ld1 dp;
    
    if (k0 == 1 && l0 == 1 && n0 == 1 && crystaltype == "si")
    { 
        dp = 3.13559;
        cpot = {
            0.00083547944310269, -0.0009866257240695535, 7.85849577178071 * pow(10, -18),
            0.002625893904708435, -0.0059639461516799285, 0.006736787164748017,
            7.444278911593285 * pow(10, -18), -0.017055968641526673, 0.03815261022457672,
            -0.042275372011778896, 2.4874088202401234 * pow(10, -16), 0.09886510785275468,
            -0.20796208250074857, 0.21683477308524746, -2.5140442684701155 * pow(10, -16),
            -0.4903539867567273, 1.0690479804477997, -1.2011945573072613,
            5.929927972277037 * pow(10, -16), 5.433460203849931, -13.90674603792947,
            5.433460203849931, 5.929927972277037 * pow(10, -16), -1.2011945573072613,
            1.0690479804477997, -0.4903539867567273, -2.5140442684701155 * pow(10, -16),
            0.21683477308524746, -0.20796208250074857, 0.09886510785275468,
            2.4874088202401234 * pow(10, -16), -0.042275372011778896, 0.03815261022457672,
            -0.017055968641526673, 7.444278911593285 * pow(10, -18), 0.006736787164748017,
            -0.0059639461516799285, 0.002625893904708435, 7.85849577178071 * pow(10, -18),
            -0.0009866257240695535, 0.00083547944310269};
        return {dp, cpot};
    }
    else if (k0 == 1 && l0 == 0 && n0 == 0 && crystaltype == "si")
    { 
        dp = 1.35775;
        cpot = {
            -1.9856702973933365 * pow(10, -13),
            2.8296527125853387 * pow(10, -12), -3.518751286601108 * pow(10, -11),
            3.81832984043055 * pow(10, -10), -3.6156614317646607 * pow(10, -9),
            2.987662901708929 * pow(10, -8), -2.1542971976402428 * pow(10, -7),
            1.355547042509427 * pow(10, -6), -7.443680319372452 * pow(10, -6),
            0.00003568419681303494, -0.00014956682386464492,
            0.0005511021652636261, -0.001813585160146616,
            0.005515945296192517, -0.01622437447835828,
            0.04683175878704749, -0.12690823070270837,
            0.31535715897428845, -0.8172178751820525,
            2.574630141524834, -13.90674603792947,
            2.574630141524834, -0.8172178751820525,
            0.31535715897428845, -0.12690823070270837,
            0.04683175878704749, -0.01622437447835828,
            0.005515945296192517, -0.001813585160146616,
            0.0005511021652636261, -0.00014956682386464492,
            0.00003568419681303494, -7.443680319372452 * pow(10, -6),
            1.355547042509427 * pow(10, -6), -2.1542971976402428 * pow(10, -7),
            2.987662901708929 * pow(10, -8), -3.6156614317646607 * pow(10, -9),
            3.81832984043055 * pow(10, -10), -3.518751286601108 * pow(10, -11),
            2.8296527125853387 * pow(10, -12), -1.9856702973933365 * pow(10, -13)};
        return {dp, cpot};
    }
    else if (k0 == 1 && l0 == 1 && n0 == 0 && crystaltype == "si")
    { 
        dp = 1.9201484643120699; 
        cpot = {
            -1.6404546940030496 * pow (10, -7), 6.192693758491359 * pow (10, -7), 
            -2.1838191384472613 * pow (10, -6), 7.194378429083061 * pow (10, -6), 
            -0.000022144503619153267, 
            0.00006370867326499635, -0.00017148442178514206, 
            0.00043289629400460457, -0.0010301629879408143, 
            0.002333234584450533, -0.005104842545865061, 
            0.010973763551363667, -0.023419159002825176, 
            0.04933841847346377, -0.10054227590506841, 
            0.19496849324298202, -0.3684800540426591, 
            0.7241983004629438, -1.4880803303113102, 
            4.339425523097554, -13.90674603792947, 
            4.339425523097554, -1.4880803303113102, 
            0.7241983004629438, -0.3684800540426591, 
            0.19496849324298202, -0.10054227590506841, 
            0.04933841847346377, -0.023419159002825176, 
            0.010973763551363667, -0.005104842545865061, 
            0.002333234584450533, -0.0010301629879408143, 
            0.00043289629400460457, -0.00017148442178514206, 
            0.00006370867326499635, -0.000022144503619153267, 
            7.194378429083061 * pow (10, -6), -2.1838191384472613 * pow (10, -6), 
            6.192693758491359 * pow (10, -7), -1.6404546940030496 * pow (10, -7)};
        return {dp, cpot};
    }
    else if (k0 == 1 && l0 == 1 && n0 == 0 && crystaltype == "c") 
    {   dp = 1.26105; 
        cpot = {
            -1.1276785549097731 * pow (10, -7), 4.921678804393794 * pow (10, -7), 
            -1.9917021674639757 * pow (10, -6), 7.47340824102231 * pow (10, -6), 
            -0.000026001368090562436, 
            0.00008387970116727148, -0.0002509000407927266, 
            0.00069586945159228, -0.0017895259230571127, 
            0.004267086198546468, -0.009434272956009272, 
            0.019340866785880345, -0.036772195689235135, 
            0.06495371898102878, -0.10777342296702845, 
            0.1757004580026742, -0.30978621881507457, 
            0.6151475757872573, -1.3733210655843027, 
            5.081355664321377, -21.172653544931986, 
            5.081355664321377, -1.3733210655843027, 
            0.6151475757872573, -0.30978621881507457, 
            0.1757004580026742, -0.10777342296702845, 
            0.06495371898102878, -0.036772195689235135, 
            0.019340866785880345, -0.009434272956009272, 
            0.004267086198546468, -0.0017895259230571127, 
            0.00069586945159228, -0.0002509000407927266, 
            0.00008387970116727148, -0.000026001368090562436, 
            7.47340824102231 * pow (10, -6), -1.9917021674639757 * pow (10, -6), 
            4.921678804393794 * pow (10, -7), -1.1276785549097731 * pow (10, -7)};
        return {dp, cpot};
    }
    else if (k0 == 1 && l0 == 1 && n0 == 1 && crystaltype == "c")
    {
        dp = 2.05929; 
        cpot = {
            0.0014265517812256272, -0.0017528531127072096, 1.436066270894863 * pow (10, -17), 
            0.004861526846768226, -0.01097343139938665, 0.012038691030198246, 
            1.2594348744032608 * pow (10, -17), -0.026623839513234278, 0.053702384403747125, 
            -0.05276596848913778, 2.7345014937951195 * pow (10, -16), 0.09648674437625458,
            -0.18507041876954172, 0.18370378581364363, -2.1132041710383178 * pow (10, -16 ),
            -0.41605249593734633, 0.9268221169602896, -1.1503407781036568, 
            6.559752649993274 * pow (10, -16), 7.150978979234087, -21.172653544931986, 
            7.150978979234087, 6.559752649993274 * pow (10, -16), -1.1503407781036568, 
            0.9268221169602896, -0.41605249593734633, -2.1132041710383178 * pow (10, -16), 
            0.18370378581364363, -0.18507041876954172, 0.09648674437625458, 
            2.7345014937951195 * pow (10, -16), -0.05276596848913778, 0.053702384403747125,
            -0.026623839513234278, 1.2594348744032608 * pow (10, -17), 0.012038691030198246, 
            -0.01097343139938665, 0.004861526846768226, 1.436066270894863 * pow (10, -17), 
            -0.0017528531127072096 , 0.0014265517812256272};
        return {dp, cpot};
    }
    else if (k0 == 1 && l0 == 0 && n0 == 0 && crystaltype == "c") 
    {
        dp = 0.8917;
        cpot = {
            -3.0819229133314926 * pow (10, -14), 5.87053476047013 * pow (10, -13), 
            -9.613907087940193 * pow (10, -12), 1.3535943778556432 * pow (10, -10), 
            -1.6384889905176582 * pow (10, -9), 1.7051589849553586 * pow (10, -8), 
            -1.525641847022871 * pow (10, -7), 1.1735637109397896 * pow (10, -6), 
            -7.761168555680368 * pow (10, -6), 
            0.000044127967404113623, -0.00021570841933776306, 
            0.0009065392946974602, -0.0032754618220711238, 
            0.010174789107963532, -0.027175380694634634, 
            0.06252375498323556, -0.12719532677318163, 
            0.2667046510796767, -0.6966493258250266, 
            2.731742703070343, -21.172653544931986, 
            2.731742703070343, -0.6966493258250266, 
            0.2667046510796767, -0.12719532677318163, 
            0.06252375498323556, -0.027175380694634634, 
            0.010174789107963532, -0.0032754618220711238, 
            0.0009065392946974602, -0.00021570841933776306, 
            0.000044127967404113623, -7.761168555680368 * pow (10, -6), 
            1.1735637109397896 * pow (10, -6), -1.525641847022871 * pow (10, -7), 
            1.7051589849553586 * pow (10, -8), -1.6384889905176582 * pow (10, -9), 
            1.3535943778556432 * pow (10, -10), -9.613907087940193 * pow (10, -12), 
            5.87053476047013 * pow (10, -13), -3.0819229133314926 * pow (10, -14)};
        return {dp, cpot};
    }
    else if (k0 == 1 && l0 == 1 && n0 == 0 && crystaltype == "ge") 
    {
        dp = 2.00041; 
        cpot = {
            -2.4542031328897046 * pow (10, -8), 1.2537373259958535 * pow (10, -7), 
            -5.890874294867771 * pow (10, -7), 2.5458349437869585 * pow (10, -6), 
            -0.000010119494672875412, 
            0.00003699727050436826, -0.00012441765037822183, 
            0.000384922604613616, -0.0010962777931523499, 
            0.002879951239617975, -0.0070158645813649285, 
            0.016040521508906818, -0.035157333519642606, 
            0.07577692235986949, -0.1625441575141206, 
            0.3414705157333183, -0.6877103073149415, 
            1.3664232474961995, -2.7775279069447194, 
            6.4225470918680925, -15.58906388110997, 
            6.4225470918680925, -2.7775279069447194, 
            1.3664232474961995, -0.6877103073149415, 
            0.3414705157333183, -0.1625441575141206, 
            0.07577692235986949, -0.035157333519642606, 
            0.016040521508906818, -0.0070158645813649285, 
            0.002879951239617975, -0.0010962777931523499, 
            0.000384922604613616, -0.00012441765037822183, 
            0.00003699727050436826, -0.000010119494672875412, 
            2.5458349437869585 * pow (10, -6), -5.890874294867771 * pow (10, -7), 
            1.2537373259958535 * pow (10, -7), -2.4542031328897046 * pow (10, -8)};
        return {dp, cpot};
    }
    else if (k0 == 1 && l0 == 1 && n0 == 1 && crystaltype == "ge") 
    {
        dp = 3.26665; 
        cpot = {
            0.0008526272706297835, -0.0011129232992799216, 
            9.671312133754144 * pow (10, -18), 0.0034722245752584163, 
            -0.008329618979087235, 0.009767128808026077, 
            1.103738148431546 * pow (10, -17), -0.02562313993385259, 
            0.05802065275418782, -0.06557790532759492, 
            3.987873311842115 * pow (10, -16), 0.16592413166778092, 
            -0.36723078753144345, 0.3994607965064238, 
            -1.998430820126219 * pow (10, -15), -0.9253159416928355, 
            2.017137403378605, -2.2089909276322612, 
            5.5161341874624214 * pow (10, -15), 6.955171953632021, 
            -15.58906388110997, 6.955171953632021, 
            5.5161341874624214 * pow (10, -15), -2.2089909276322612, 
            2.017137403378605, -0.9253159416928355, 
            -1.998430820126219 * pow (10, -15), 0.3994607965064238,
            -0.36723078753144345, 0.16592413166778092, 
            3.987873311842115 * pow (10, -16), -0.06557790532759492, 
            0.05802065275418782, -0.02562313993385259, 
            1.103738148431546 * pow (10, -17), 0.009767128808026077, 
            -0.008329618979087235, 0.0034722245752584163, 
            9.671312133754144 * pow (10, -18), -0.0011129232992799216, 
            0.0008526272706297835};
        return {dp, cpot};
    }

    else if (k0 == 1 && l0 == 0 && n0 == 0 && crystaltype == "ge") 
    {
        dp = 1.4145; 
        cpot = {
            -1.3345525371436056 * pow (10, -15), 3.48279383638218 * pow (10, -14),
            -7.689075243520236 * pow (10, -13), 1.4360670121597631 * pow (10, -11),
            -2.2689741766336658 * pow (10, -10) , 3.0327642493123866 * pow (10, -9), 
            -3.42927254932831 * pow (10, -8), 3.2803457361144563 * pow (10, -7), 
            -2.654554786197116 * pow (10, -6),
            0.00001817274628978669, -0.0001052515818142472, 
            0.0005159224261672529, -0.0021452328355364217, 
            0.007643665794082016, -0.024110331383628877, 
            0.07176267589339667, -0.21050226010200654, 
            0.58231592107331, -1.5416057395581846, 
            4.344839622082629, -15.58906388110997, 
            4.344839622082629, -1.5416057395581846, 
            0.58231592107331, -0.21050226010200654, 
            0.07176267589339667, -0.024110331383628877, 
            0.007643665794082016, -0.0021452328355364217, 
            0.0005159224261672529, -0.0001052515818142472, 
            0.00001817274628978669, -2.654554786197116 * pow (10, -6), 
            3.2803457361144563 * pow (10, -7), -3.42927254932831 * pow (10, -8), 
            3.0327642493123866 * pow (10, -9), -2.2689741766336658 * pow (10, -10), 
            1.4360670121597631 * pow (10, -11), -7.689075243520236 * pow (10, -13), 
            3.48279383638218 * pow (10, -14), -1.3345525371436056 * pow (10, -15)};
        return {dp, cpot};
    }
    else if (k0 == 1 && l0 == 1 && n0 == 0 && crystaltype == "w") 
    {
        dp = 2.3334523779156067; 
        cpot = {
            -0.014550833927849083, -0.021926004705406607, -0.032354671565973984, 
            -0.046761277622265306, -0.06621747478371258, -0.091955420484176,
            -0.12546204066266725, -0.16879100349942325, -0.22532138673354637,
            -0.3012322907124337, -0.4077881899634045, -0.5639289964279592, 
            -0.7975956836483334, -1.143471807620767, -1.6384541420542922,
            -2.3381118703258825, -3.4204250001995535, -5.351101288430016, 
            -8.916001272180678, -17.253360850252008, -33.35588375745145, 
            -17.253360850252008, -8.916001272180678, -5.351101288430016, 
            -3.4204250001995535, -2.3381118703258825, -1.6384541420542922, 
            -1.143471807620767, -0.7975956836483334, -0.5639289964279592, 
            -0.4077881899634045, -0.3012322907124337, -0.22532138673354637, 
            -0.16879100349942325, -0.12546204066266725, -0.091955420484176,
            -0.06621747478371258, -0.046761277622265306, -0.032354671565973984, 
            -0.021926004705406607,-0.014550833927849083};
        return {dp, cpot};
    }
    else if (k0 == 1 && l0 == 1 && n0 == 1 && crystaltype == "w") 
    {
        dp = 0.9526279441628825; 
        cpot = {
            -1.076245582650239 * pow (10, -11), -1.2596061056095884 * pow (10, -10), 
            -1.299488133862719 * pow (10, -9), -1.1817457934884241 * pow (10, -8), 
            -9.47305367153397 * pow (10, -8), -6.693759686000083 * pow (10, -7), 
            -4.169312245098132 * pow (10, -6), -0.000022891429384758735, 
            -0.00011078873105041371, -0.0004726421515757212, -0.001777393195846447, 
            -0.005891816966833513, -0.01721632565605005, -0.04436339322583067, 
            -0.10121085753470994, -0.20983402765244152, -0.43456554030122346, 
            -1.0079292205659602, -2.4249832314909274, -7.024614278484769, 
            -33.35588375745145, -7.024614278484769, -2.4249832314909274, 
            -1.0079292205659602, -0.43456554030122346, -0.20983402765244152, 
            -0.1012108575347099, -0.04436339322583067, -0.01721632565605005, 
            -0.005891816966833513, -0.001777393195846447, -0.0004726421515757212, 
            -0.00011078873105041371, -0.000022891429384758735, -4.169312245098132 * pow (10, -6),
            -6.693759686000083 * pow (10, -7), -9.47305367153397 * pow (10, -8), 
            -1.1817457934884241 * pow (10, -8), -1.299488133862719 * pow (10, -9), 
            -1.2596061056095884 * pow (10, -10), -1.076245582650239 * pow (10, -11)};
        return {dp, cpot};
    }
    else if (k0 == 1 && l0 == 0 && n0 == 0 && crystaltype == "w") 
    {
        dp = 1.65; 
        cpot = {
            -0.00021711338412334374, -0.0004929403722600468, -0.0010731001248195282,
            -0.002239877086450769, -0.004482769043168181, -0.008602164677756236,
            -0.015827556666879503, -0.027924925747252482, -0.04725626060193743, 
            -0.07678050131633864, -0.12015785136047805, -0.1826717142298991,
            -0.2747802793141709, -0.42084724659162154, -0.6724307161564471, 
            -1.1144146916354174, -1.8508045331939567, -3.1030337067794043, 
            -5.815454027478652, -12.792748125972803, -33.35588375745145, 
            -12.792748125972803, -5.815454027478652, -3.1030337067794043,
            -1.8508045331939567, -1.1144146916354174, -0.6724307161564471, 
            -0.42084724659162154, -0.2747802793141709, -0.1826717142298991, 
            -0.12015785136047805, -0.07678050131633864, -0.04725626060193743, 
            -0.027924925747252482, -0.015827556666879503, -0.008602164677756236, 
            -0.004482769043168181, -0.002239877086450769, -0.0010731001248195282, 
            -0.0004929403722600468, -0.00021711338412334374};
        return {dp, cpot};
    }
	else 
	{
		dp = 3.13559;
        cpot = {
            0.00083547944310269, -0.0009866257240695535, 7.85849577178071 * pow(10, -18),
            0.002625893904708435, -0.0059639461516799285, 0.006736787164748017,
            7.444278911593285 * pow(10, -18), -0.017055968641526673, 0.03815261022457672,
            -0.042275372011778896, 2.4874088202401234 * pow(10, -16), 0.09886510785275468,
            -0.20796208250074857, 0.21683477308524746, -2.5140442684701155 * pow(10, -16),
            -0.4903539867567273, 1.0690479804477997, -1.2011945573072613,
            5.929927972277037 * pow(10, -16), 5.433460203849931, -13.90674603792947,
            5.433460203849931, 5.929927972277037 * pow(10, -16), -1.2011945573072613,
            1.0690479804477997, -0.4903539867567273, -2.5140442684701155 * pow(10, -16),
            0.21683477308524746, -0.20796208250074857, 0.09886510785275468,
            2.4874088202401234 * pow(10, -16), -0.042275372011778896, 0.03815261022457672,
            -0.017055968641526673, 7.444278911593285 * pow(10, -18), 0.006736787164748017,
            -0.0059639461516799285, 0.002625893904708435, 7.85849577178071 * pow(10, -18),
            -0.0009866257240695535, 0.00083547944310269};
        return {dp, cpot};
	}
}