## Code to generate a DCN cell.

class DCN(object):

    def __init__(self,fiberD = 2.0,**kwargs):
        from neuron import h
        self.fiberD = fiberD
        self.model_globals()
        self.initcell()


# //Variables//
    def model_globals(self):
        import math
        from neuron import h
        h.celsius = 36    #//30 degC
        h.v_init = -60    #//mV

        # //topological parameters
        self.dendelements = 60
        self.somaelements = 1
        self.iselements = 1       #//initial segment elements
        self.axonnodes = 27     
        self.paranodes1 = 52
        self.paranodes2 = 52
        self.axoninter = 78
        # self.total = dendelements+somaelements+iselements+axonnodes+paranodes1+paranodes2+axoninter

        # //morphological parameters of the axon
        # fiberD=2            //fiber diameter
        self.nodelength = 1      #//node length
        self.paralength1 = 3     #//MYSA length
        self.paralength2 = 10    #//FLUT length
        self.space_p1 = 0.002  
        self.space_p2 = 0.004
        self.space_i = 0.004

        # //electrical parameters of the axon
        self.rhoa = 0.7e6     #// Ohm-um cytoplasmic resistance
        self.mycm = 0.1       #// uF/cm2 lamella membrane capacitance
        self.mygm = 0.001     #// S/cm2 lamella membrane conductance
    
        if self.fiberD==1.0:
            self.axonD=0.8
            self.nodeD=0.7
            self.paraD1=0.7
            self.paraD2=0.8
            self.deltax=200.0
            self.nl=20

        if self.fiberD==2.0:
            self.axonD=1.6
            self.nodeD=1.4
            self.paraD1=1.4
            self.paraD2=1.6 
            self.deltax=200.0
            self.nl=30

        if self.fiberD == 5.7:
            self.g=0.605 
            self.axonD=3.4 
            self.nodeD=1.9
            self.paraD1=1.9 
            self.paraD2=3.4 
            self.deltax=500.
            self.paralength2=35.
            self.nl=80.

        self.Rpn0 = (self.rhoa*.01)/(math.pi*((((self.nodeD/2.)+self.space_p1)**2)-((self.nodeD/2.)**2)))
        self.Rpn1 = (self.rhoa*.01)/(math.pi*((((self.paraD1/2.)+self.space_p1)**2)-((self.paraD1/2.)**2)))
        self.Rpn2 = (self.rhoa*.01)/(math.pi*((((self.paraD2/2.)+self.space_p2)**2)-((self.paraD2/2.)**2)))
        self.Rpx  = (self.rhoa*.01)/(math.pi*((((self.axonD/2.)+self.space_i)**2)-((self.axonD/2.)**2)))

        self.interlength = (self.deltax-self.nodelength-(2*self.paralength1)-(2*self.paralength2))/6.

        self.cai0_cacum = 1e-4
        self.cai0_ca_ion = 1e-4
        self.ki0_k_ion = 140
        self.nai0_na_ion = 15 
        self.cao0_ca_ion = 2.4
        self.ko0_k_ion = 6.24
        self.nao0_na_ion = 140 
        self.my_gna_na = 0.028          
        self.my_gna_nal = 5.e-4         
        self.my_vtraub_iNaP = -65.
        self.my_gh_ih =  .0002 
        self.my_gk_kdrf = .006 
        self.my_ba = 0.001 
        self.my_gk_ikca = 0.000025
        self.my_gcal_hva = .0001 
        self.my_gcan_hva = .001 


        self.my_g_pas = 2.0e-4
        self.my_e_pas = -52.0


        self.my_g_CaT = .0001



    def initcell(self):
        from neuron import h
        self.soma = []
        for i in range(self.somaelements): 
            self.soma.append(h.Section())

            self.soma[i].nseg = 1
            self.soma[i].Ra = 150.0   
            self.soma[i].cm = 2.


            self.soma[i].insert('myions')   #//dummy mechanism to set up ion concentrations for na and k        

            self.soma[i].insert('Na')                    #//MODEL FOR NA+ CURRENT
            self.soma[i].gna_Na = self.my_gna_na      
            self.soma[i].insert('iNaP')                     #/MODEL FOR NA+ LEAK CURRENT
            self.soma[i].gnabar_iNaP = self.my_gna_nal
            self.soma[i].vtraub_iNaP = self.my_vtraub_iNaP
            self.soma[i].insert('Ih')                    #//MODEL FOR K+ INWARD RECTIFIER CURRENT
            self.soma[i].gh_Ih = self.my_gh_ih
            self.soma[i].insert('KDRf')                    #//MODEL FOR K+ DELAYED RECTIFIER CURRENT (fast-deactivating)
            self.soma[i].gk_KDRf = self.my_gk_kdrf

            self.soma[i].insert('B_A')
            self.soma[i].gkbar_B_A = self.my_ba


            self.soma[i].insert('iKCa')                    #//MODEL FOR SMALL COND CA2+ ACTIVATED K+ CURRENT
            self.soma[i].gbar_iKCa = self.my_gk_ikca
            self.soma[i].insert('cacum')                    #//MODEL FOR CA2+ ACCUMULATION
            self.soma[i].cai0_cacum = self.cai0_ca_ion
            self.soma[i].insert('HVA')                    #//MODEL FOR HVA Ca2+ CURRENT
            self.soma[i].gcaL_HVA = self.my_gcal_hva
            self.soma[i].gcaN_HVA = self.my_gcan_hva
            self.soma[i].insert('extracellular')

            self.soma[i].insert('pas')
            self.soma[i].e_pas = self.my_e_pas
            self.soma[i].g_pas = self.my_g_pas

            self.soma[i].insert('CaT')
            self.soma[i].gcaT_CaT = self.my_g_CaT

        self.dend = []
        for i in range(self.dendelements):
            self.dend.append(h.Section())
            self.dend[i].nseg = 1
            self.dend[i].Ra = 150.0 
            self.dend[i].cm = 2.

            self.dend[i].insert('myions')   #//dummy mechanism to set up ion concentrations for na and k        
            self.dend[i].insert('Na')       #                //MODEL FOR NA+ CURRENT
            self.dend[i].gna_Na = 1e-7               
            self.dend[i].insert('iNaP')     #               //MODEL FOR NA+ LEAK CURRENT
            self.dend[i].gnabar_iNaP = self.my_gna_nal
            self.dend[i].vtraub_iNaP = self.my_vtraub_iNaP
            self.dend[i].insert('Ih')     #              //MODEL FOR K+ INWARD RECTIFIER CURRENT
            self.dend[i].gh_Ih = self.my_gh_ih
            self.dend[i].insert('KDRf')  #                 //MODEL FOR FAST-DEACTIVATING K+ CURRENT (delayed rectifier)
            self.dend[i].gk_KDRf = self.my_gk_kdrf

            self.dend[i].insert('B_A')
            self.dend[i].gkbar_B_A = self.my_ba
            self.dend[i].insert('iKCa')   #                //MODEL FOR SMALL COND CA2+ ACTIVATED K+ CURRENT
            self.dend[i].gbar_iKCa = self.my_gk_ikca
            self.dend[i].insert('cacum') #                   //MODEL FOR CA2+ ACCUMULATION
            self.dend[i].cai0_cacum = self.cai0_ca_ion
            self.dend[i].insert('HVA') #                 //MODEL FOR HVA Ca2+ CURRENT
            self.dend[i].gcaN_HVA = self.my_gcan_hva
            self.dend[i].gcaL_HVA = self.my_gcal_hva
            self.dend[i].insert('extracellular')

            self.dend[i].insert('pas')
            self.dend[i].e_pas = self.my_e_pas
            self.dend[i].g_pas = self.my_g_pas

            self.dend[i].insert('CaT')
            self.dend[i].gcaT_CaT = self.my_g_CaT

            #xraxial=1e+09 xg=1e+09 xc=0
            

        self.initseg = h.Section()
        self.initseg.nseg = 1
        self.initseg.Ra = 150.0
        self.initseg.cm = 2.

        self.initseg.insert('myions')
        self.initseg.insert('Na')
        self.initseg.gna_Na = self.my_gna_na
        self.initseg.insert('iNaP')
        self.initseg.gnabar_iNaP = self.my_gna_nal
        self.initseg.vtraub_iNaP = self.my_vtraub_iNaP
        self.initseg.insert('Ih')
        self.initseg.gh_Ih = self.my_gh_ih
        self.initseg.insert('KDRf')
        self.initseg.gk_KDRf = self.my_gk_kdrf

        self.initseg.insert('B_A')
        self.initseg.gkbar_B_A = self.my_ba
        self.initseg.insert('iKCa')
        self.initseg.gbar_iKCa = self.my_gk_ikca
        self.initseg.insert('cacum')
        self.initseg.cai0_cacum = self.cai0_ca_ion
        self.initseg.insert('HVA')
        self.initseg.gcaL_HVA = self.my_gcal_hva
        self.initseg.gcaN_HVA = self.my_gcan_hva
        self.initseg.insert('extracellular')

        self.initseg.insert('pas')
        self.initseg.e_pas = self.my_e_pas
        self.initseg.g_pas = self.my_g_pas

        self.initseg.insert('CaT')
        self.dend[i].gcaT_CaT = self.my_g_CaT

        self.node = []
        for i in range(self.axonnodes - 1):
            self.node.append(h.Section())
            self.node[i].nseg = 1
            self.node[i].Ra = self.rhoa/10000.
            self.node[i].cm = 2.
            self.node[i].insert('axnode75')
            self.node[i].gnabar_axnode75 = 2.0  
            self.node[i].gnapbar_axnode75 = 0.05 
            self.node[i].gkbar_axnode75 = 0.07
            self.node[i].gl_axnode75 = 0.005
            self.node[i].ek_axnode75 = -85.
            self.node[i].ena_axnode75 = 55.
            self.node[i].el_axnode75 = -60.
            self.node[i].vshift_axnode75 = 15.
            self.node[i].vtraub_axnode75 = -80.
            self.node[i].insert('extracellular')
            self.node[i].xraxial[0]=self.Rpn0 
            self.node[i].xg[0]=1e10 
            self.node[i].xc[0]=0
        
        self.node.append(h.Section())
        self.node[-1].nseg = 1
        self.node[-1].Ra = self.rhoa/10000.0
        self.node[-1].cm = 2.0
        self.node[-1].insert('pas')
        self.node[-1].g_pas = 0.0001
        self.node[-1].e_pas = h.v_init 
        self.node[-1].insert('extracellular') 
        self.node[-1].xraxial[0]=self.Rpn0 
        self.node[-1].xg[0]=1e10 
        self.node[-1].xc[0]=0

        self.MYSA = []
        for i in range(self.paranodes1):
            self.MYSA.append(h.Section())
            self.MYSA[i].nseg = 1
            self.MYSA[i].Ra = self.rhoa/10000.0
            self.MYSA[i].cm = 2.
            self.MYSA[i].insert('pas')
            self.MYSA[i].g_pas = 0.0001 
            self.MYSA[i].e_pas = h.v_init
            self.MYSA[i].insert('extracellular')
            self.MYSA[i].xraxial[0]=self.Rpn1
            self.MYSA[i].xg[0]=self.mygm/(self.nl*2.) 
            self.MYSA[i].xc[0]=self.mycm/(self.nl*2.)
            
        self.FLUT = []
        for i in range(self.paranodes2):
            self.FLUT.append(h.Section())
            self.FLUT[i].nseg = 1
            self.FLUT[i].Ra = self.rhoa/10000.
            self.FLUT[i].cm = 2.
            self.FLUT[i].insert('parak75')
            self.FLUT[i].gkbar_parak75 = 0.02
            self.FLUT[i].ek_parak75 = -85.
            self.FLUT[i].vshift_parak75 = 15.
            self.FLUT[i].insert('pas')
            self.FLUT[i].g_pas = 0.0001     
            self.FLUT[i].e_pas = h.v_init
            self.FLUT[i].insert('extracellular')
            self.FLUT[i].xraxial[0]=self.Rpn2
            self.FLUT[i].xg[0]=self.mygm/(self.nl*2.)
            self.FLUT[i].xc[0]=self.mycm/(self.nl*2.)
        
        self.STIN = []
        for i in range(self.axoninter):
            self.STIN.append(h.Section())
            self.STIN[i].nseg = 1
            self.STIN[i].Ra = self.rhoa/10000.
            self.STIN[i].cm = 2.
            self.STIN[i].insert('pas')
            self.STIN[i].g_pas = 0.0001
            self.STIN[i].e_pas = h.v_init
            self.STIN[i].insert('extracellular')
            self.STIN[i].xraxial[0]=self.Rpx
            self.STIN[i].xg[0]=self.mygm/(self.nl*2.) 
            self.STIN[i].xc[0]=self.mycm/(self.nl*2.)
            


        for i in range(self.somaelements - 1):
            self.soma[i+1].connect(self.soma[i](1))

        self.dend[0].connect(self.soma[0](1))
        self.dend[1].connect(self.dend[0](1))
        self.dend[2].connect(self.dend[1](1))
        self.dend[3].connect(self.dend[2](1))
        self.dend[4].connect(self.dend[3](1))
        self.dend[5].connect(self.dend[4](1))
        self.dend[6].connect(self.dend[5](1))
        self.dend[7].connect(self.dend[6](1))
        self.dend[8].connect(self.dend[7](1))
        self.dend[9].connect(self.dend[8](1))
        self.dend[10].connect(self.dend[9](1))
        self.dend[11].connect(self.dend[10](1))
        self.dend[12].connect(self.dend[11](1))
        self.dend[13].connect(self.dend[12](1))
        self.dend[14].connect(self.dend[13](1))
        self.dend[15].connect(self.dend[6](1))    #//branch
        self.dend[16].connect(self.dend[15](1))
        self.dend[17].connect(self.dend[0](1))    #//branch
        self.dend[18].connect(self.dend[17](1))
        self.dend[19].connect(self.dend[18](1))
        self.dend[20].connect(self.dend[19](1))
        self.dend[21].connect(self.dend[20](1))
        self.dend[22].connect(self.dend[21](1))
        self.dend[23].connect(self.dend[22](1))
        self.dend[24].connect(self.dend[23](1))
        self.dend[25].connect(self.dend[0](1))    #//branch
        self.dend[26].connect(self.dend[25](1))
        self.dend[27].connect(self.dend[26](1))
        self.dend[28].connect(self.dend[27](1))
        self.dend[29].connect(self.dend[28](1))
        self.dend[30].connect(self.dend[29](1))
        self.dend[31].connect(self.dend[30](1))
        self.dend[32].connect(self.dend[31](1))
        self.dend[33].connect(self.dend[32](1))
        self.dend[34].connect(self.dend[33](1))
        self.dend[35].connect(self.dend[34](1))
        self.dend[36].connect(self.dend[35](1))
        self.dend[37].connect(self.dend[36](1))
        self.dend[38].connect(self.dend[28](1))    #//branch
        self.dend[39].connect(self.dend[38](1))
        self.dend[40].connect(self.dend[39](1))
        self.dend[41].connect(self.soma[0](0))    #//branch
        self.dend[42].connect(self.dend[41](1))
        self.dend[43].connect(self.dend[42](1))
        self.dend[44].connect(self.dend[43](1))
        self.dend[45].connect(self.dend[44](1))
        self.dend[46].connect(self.dend[45](1))
        self.dend[47].connect(self.dend[46](1))
        self.dend[48].connect(self.dend[47](1))
        self.dend[49].connect(self.dend[48](1))
        self.dend[50].connect(self.dend[49](1))
        self.dend[51].connect(self.dend[50](1))
        self.dend[52].connect(self.dend[41](1))    #//branch
        self.dend[53].connect(self.dend[52](1))
        self.dend[54].connect(self.dend[53](1))
        self.dend[55].connect(self.dend[54](1))
        self.dend[56].connect(self.dend[55](1))
        self.dend[57].connect(self.dend[56](1))
        self.dend[58].connect(self.dend[47](1))
        self.dend[59].connect(self.dend[58](1))
        self.initseg.connect(self.soma[0](1))
        self.node[0].connect(self.initseg(1))

        for i in range(self.axonnodes-1):
            self.MYSA[2*i].connect(self.node[i](1))
            self.FLUT[2*i].connect(self.MYSA[2*i](1))
            self.STIN[3*i].connect(self.FLUT[2*i](1))
            self.STIN[3*i+1].connect(self.STIN[3*i](1))
            self.STIN[3*i+2].connect(self.STIN[3*i+1](1))
            self.FLUT[2*i+1].connect(self.STIN[3*i+2](1))
            self.MYSA[2*i+1].connect(self.FLUT[2*i+1](1))
            self.node[i+1].connect(self.MYSA[2*i+1](1))


        # 
        h.pt3dadd(-8515.159,5268.953,-1303.344,4, sec = self.soma[0])
        h.pt3dadd(-8511.8016,5267.0221,-1306.5062,10, sec = self.soma[0])

        h.pt3dadd(-8508.4442,5265.0913,-1309.6684,18, sec = self.soma[0])
        h.pt3dadd(-8505.0868,5263.1604,-1312.8306,18, sec = self.soma[0])
        h.pt3dadd(-8501.7294,5261.2295,-1315.9929,18, sec = self.soma[0])
        h.pt3dadd(-8498.7863,5262.1054,-1319.9388,10, sec = self.soma[0])
        h.pt3dadd(-8495.8431,5262.9813,-1323.8848,4, sec = self.soma[0])
        
        h.pt3dadd(-8494.9405,5271.2959,-1329.367,4.0, sec = self.dend[0])
        h.pt3dadd(-8484.0863,5275.1863,-1340.6829,1.53, sec = self.dend[0])

        h.pt3dadd(-8484.0863,5275.1863,-1340.6829,1.53, sec = self.dend[1])
        h.pt3dadd(-8476.5526,5279.9636,-1345.2019,1.53, sec = self.dend[1])

        h.pt3dadd(-8476.5526,5279.9636,-1345.2019,1.53, sec = self.dend[2])
        h.pt3dadd(-8469.0189,5284.7409,-1349.7208,1.53, sec = self.dend[2])

        h.pt3dadd(-8469.0189,5284.7409,-1349.7208,1.53, sec = self.dend[3])
        h.pt3dadd(-8461.4851,5289.5182,-1354.2397,1.53, sec = self.dend[3])

        h.pt3dadd(-8461.4851,5289.5182,-1354.2397,1.53, sec = self.dend[4])
        h.pt3dadd(-8458.6422,5298.1982,-1358.311,1.53, sec = self.dend[4])

        h.pt3dadd(-8458.6422,5298.1982,-1358.311,1.53, sec = self.dend[5])
        h.pt3dadd(-8455.7993,5306.8782,-1362.3823,1.53, sec = self.dend[5])

        h.pt3dadd(-8455.7993,5306.8782,-1362.3823,1.53, sec = self.dend[6])
        h.pt3dadd(-8451.4257,5313.8504,-1368.062,1.53, sec = self.dend[6])

        h.pt3dadd(-8451.4257,5313.8504,-1368.062,1.53, sec = self.dend[7])
        h.pt3dadd(-8444.985,5318.6494,-1374.0191,1.53, sec = self.dend[7])

        h.pt3dadd(-8444.985,5318.6494,-1374.0191,1.53, sec = self.dend[8])
        h.pt3dadd(-8438.5444,5323.4485,-1379.9762,1.53, sec = self.dend[8])

        h.pt3dadd(-8438.5444,5323.4485,-1379.9762,1.53, sec = self.dend[9])
        h.pt3dadd(-8432.1037,5328.2475,-1385.9334,1.53, sec = self.dend[9])

        h.pt3dadd(-8432.1037,5328.2475,-1385.9334,1.53, sec = self.dend[10])
        h.pt3dadd(-8428.1464,5336.1494,-1390.6131,1.53, sec = self.dend[10])

        h.pt3dadd(-8428.1464,5336.1494,-1390.6131,1.53, sec = self.dend[11])
        h.pt3dadd(-8424.189,5344.0513,-1395.2928,1.53, sec = self.dend[11])

        h.pt3dadd(-8424.189,5344.0513,-1395.2928,1.53, sec = self.dend[12])
        h.pt3dadd(-8420.2316,5351.9531,-1399.9725,1.53, sec = self.dend[12])

        h.pt3dadd(-8420.2316,5351.9531,-1399.9725,1.53, sec = self.dend[13])
        h.pt3dadd(-8413.2853,5357.368,-1404.7083,1.53, sec = self.dend[13])

        h.pt3dadd(-8413.2853,5357.368,-1404.7083,1.53, sec = self.dend[14])
        h.pt3dadd(-8406.339,5362.783,-1409.444,1.53, sec = self.dend[14])

        h.pt3dadd(-8451.4257,5313.8504,-1368.062,1.53, sec = self.dend[15])
        h.pt3dadd(-8458.4573,5316.8861,-1361.6324,1.53, sec = self.dend[15])

        h.pt3dadd(-8458.4573,5316.8861,-1361.6324,1.53, sec = self.dend[16])
        h.pt3dadd(-8465.489,5319.9219,-1355.2028,1.53, sec = self.dend[16])

        h.pt3dadd(-8484.0863,5275.1863,-1340.6829,1.53, sec = self.dend[17])
        h.pt3dadd(-8475.3193,5276.9344,-1345.1644,1.53, sec = self.dend[17])

        h.pt3dadd(-8475.3193,5276.9344,-1345.1644,1.53, sec = self.dend[18])
        h.pt3dadd(-8466.5523,5278.6826,-1349.6458,1.02, sec = self.dend[18])

        h.pt3dadd(-8466.5523,5278.6826,-1349.6458,1.02, sec = self.dend[19])
        h.pt3dadd(-8457.7853,5280.4307,-1354.1273,1.02, sec = self.dend[19])

        h.pt3dadd(-8457.7853,5280.4307,-1354.1273,1.02, sec = self.dend[20])
        h.pt3dadd(-8449.0182,5282.1789,-1358.6088,1.02, sec = self.dend[20])

        h.pt3dadd(-8449.0182,5282.1789,-1358.6088,1.02, sec = self.dend[21])
        h.pt3dadd(-8440.2512,5283.9271,-1363.0902,1.02, sec = self.dend[21])

        h.pt3dadd(-8440.2512,5283.9271,-1363.0902,1.02, sec = self.dend[22])
        h.pt3dadd(-8432.6846,5279.4051,-1367.8123,1.02, sec = self.dend[22])

        h.pt3dadd(-8432.6846,5279.4051,-1367.8123,1.02, sec = self.dend[23])
        h.pt3dadd(-8425.118,5274.8831,-1372.5344,1.02, sec = self.dend[23])

        h.pt3dadd(-8425.118,5274.8831,-1372.5344,1.02, sec = self.dend[24])
        h.pt3dadd(-8417.5514,5270.3611,-1377.2565,1.02, sec = self.dend[24])

        h.pt3dadd(-8484.0863,5275.1863,-1340.6829,1.53, sec = self.dend[25])
        h.pt3dadd(-8480.4223,5267.9191,-1346.4935,1.53, sec = self.dend[25])

        h.pt3dadd(-8480.4223,5267.9191,-1346.4935,1.53, sec = self.dend[26])
        h.pt3dadd(-8476.7584,5260.6518,-1352.3041,1.53, sec = self.dend[26])

        h.pt3dadd(-8476.7584,5260.6518,-1352.3041,1.53, sec = self.dend[27])
        h.pt3dadd(-8471.18,5253.9316,-1357.1745,1.53, sec = self.dend[27])

        h.pt3dadd(-8471.18,5253.9316,-1357.1745,1.53, sec = self.dend[28])
        h.pt3dadd(-8464.0987,5247.8658,-1360.7884,1.53, sec = self.dend[28])

        h.pt3dadd(-8464.0987,5247.8658,-1360.7884,1.53, sec = self.dend[29])
        h.pt3dadd(-8457.0173,5241.8,-1364.4023,1.53, sec = self.dend[29])

        h.pt3dadd(-8457.0173,5241.8,-1364.4023,1.53, sec = self.dend[30])
        h.pt3dadd(-8449.9359,5235.7342,-1368.0162,1.53, sec = self.dend[30])

        h.pt3dadd(-8449.9359,5235.7342,-1368.0162,1.53, sec = self.dend[31])
        h.pt3dadd(-8442.4074,5229.4631,-1366.0177,1.53, sec = self.dend[31])

        h.pt3dadd(-8442.4074,5229.4631,-1366.0177,1.53, sec = self.dend[32])
        h.pt3dadd(-8434.879,5223.1919,-1364.0191,1.53, sec = self.dend[32])

        h.pt3dadd(-8434.879,5223.1919,-1364.0191,1.53, sec = self.dend[33])
        h.pt3dadd(-8427.3505,5216.9207,-1362.0205,1.53, sec = self.dend[33])

        h.pt3dadd(-8427.3505,5216.9207,-1362.0205,1.53, sec = self.dend[34])
        h.pt3dadd(-8417.5681,5216.2068,-1360.0723,1.53, sec = self.dend[34])

        h.pt3dadd(-8417.5681,5216.2068,-1360.0723,1.53, sec = self.dend[35])
        h.pt3dadd(-8410.8479,5208.8038,-1359.8867,1.53, sec = self.dend[35])

        h.pt3dadd(-8410.8479,5208.8038,-1359.8867,1.53, sec = self.dend[36])
        h.pt3dadd(-8404.1277,5201.4009,-1359.7011,1.53, sec = self.dend[36])

        h.pt3dadd(-8404.1277,5201.4009,-1359.7011,1.53, sec = self.dend[37])
        h.pt3dadd(-8397.4074,5193.9979,-1359.5155,1.53, sec = self.dend[37])

        h.pt3dadd(-8464.0987,5247.8658,-1360.7884,1.53, sec = self.dend[38])
        h.pt3dadd(-8462.1739,5248.7294,-1370.5633,1.53, sec = self.dend[38])

        h.pt3dadd(-8462.1739,5248.7294,-1370.5633,1.53, sec = self.dend[39])
        h.pt3dadd(-8460.249,5249.5931,-1380.3383,1.02, sec = self.dend[39])

        h.pt3dadd(-8460.249,5249.5931,-1380.3383,1.02, sec = self.dend[40])
        h.pt3dadd(-8455.5809,5252.2443,-1388.7751,1.02, sec = self.dend[40])

        h.pt3dadd(-8511.8016,5267.0221,-1306.5062,4.0, sec = self.dend[41])
        h.pt3dadd(-8525.5135,5270.3276,-1298.5091,0.605, sec = self.dend[41])

        h.pt3dadd(-8525.5135,5270.3276,-1298.5091,0.605, sec = self.dend[42])
        h.pt3dadd(-8531.0271,5266.6551,-1291.0181,0.605, sec = self.dend[42])

        h.pt3dadd(-8531.0271,5266.6551,-1291.0181,0.605, sec = self.dend[43])
        h.pt3dadd(-8536.5406,5262.9827,-1283.5272,0.605, sec = self.dend[43])

        h.pt3dadd(-8536.5406,5262.9827,-1283.5272,0.605, sec = self.dend[44])
        h.pt3dadd(-8542.0541,5259.3102,-1276.0363,0.605, sec = self.dend[44])

        h.pt3dadd(-8542.0541,5259.3102,-1276.0363,0.605, sec = self.dend[45])
        h.pt3dadd(-8547.5677,5255.6378,-1268.5453,0.605, sec = self.dend[45])

        h.pt3dadd(-8547.5677,5255.6378,-1268.5453,0.605, sec = self.dend[46])
        h.pt3dadd(-8553.0812,5251.9653,-1261.0544,0.605, sec = self.dend[46])

        h.pt3dadd(-8553.0812,5251.9653,-1261.0544,0.605, sec = self.dend[47])
        h.pt3dadd(-8559.1384,5250.9453,-1253.1633,0.605, sec = self.dend[47])

        h.pt3dadd(-8559.1384,5250.9453,-1253.1633,0.605, sec = self.dend[48])
        h.pt3dadd(-8565.1956,5249.9253,-1245.2721,0.605, sec = self.dend[48])

        h.pt3dadd(-8565.1956,5249.9253,-1245.2721,0.605, sec = self.dend[49])
        h.pt3dadd(-8571.2528,5248.9053,-1237.381,0.605, sec = self.dend[49])

        h.pt3dadd(-8571.2528,5248.9053,-1237.381,0.605, sec = self.dend[50])
        h.pt3dadd(-8577.5846,5249.0893,-1229.6431,0.605, sec = self.dend[50])

        h.pt3dadd(-8577.5846,5249.0893,-1229.6431,0.605, sec = self.dend[51])
        h.pt3dadd(-8583.9164,5249.2732,-1221.9053,0.605, sec = self.dend[51])

        h.pt3dadd(-8525.5135,5270.3276,-1298.5091,0.605, sec = self.dend[52])
        h.pt3dadd(-8527.7642,5262.888,-1292.2172,0.605, sec = self.dend[52])

        h.pt3dadd(-8527.7642,5262.888,-1292.2172,0.605, sec = self.dend[53])
        h.pt3dadd(-8530.0148,5255.4484,-1285.9253,0.605, sec = self.dend[53])

        h.pt3dadd(-8530.0148,5255.4484,-1285.9253,0.605, sec = self.dend[54])
        h.pt3dadd(-8532.2655,5248.0089,-1279.6335,0.605, sec = self.dend[54])

        h.pt3dadd(-8532.2655,5248.0089,-1279.6335,0.605, sec = self.dend[55])
        h.pt3dadd(-8534.5161,5240.5693,-1273.3416,0.605, sec = self.dend[55])

        h.pt3dadd(-8534.5161,5240.5693,-1273.3416,0.605, sec = self.dend[56])
        h.pt3dadd(-8536.7668,5233.1298,-1267.0497,0.605, sec = self.dend[56])

        h.pt3dadd(-8536.7668,5233.1298,-1267.0497,0.605, sec = self.dend[57])
        h.pt3dadd(-8539.0174,5225.6902,-1260.7578,0.605, sec = self.dend[57])

        h.pt3dadd(-8539.0174,5225.6902,-1260.7578,0.605, sec = self.dend[58])
        h.pt3dadd(-8541.268,5218.2507,-1254.466,0.605, sec = self.dend[58])

        h.pt3dadd(-8541.268,5218.2507,-1254.466,0.605, sec = self.dend[59])
        h.pt3dadd(-8541.5594,5211.623,-1246.9833,0.605, sec = self.dend[59])

        h.pt3dadd(-8508.4442,5265.0913,-1309.6684,4.0, sec = self.initseg)
        h.pt3dadd(-8510.692,5263.453,-1306.016,self.nodeD, sec = self.initseg)

        h.pt3dadd(-8510.692,5263.453,-1306.016,self.nodeD, sec = self.node[0])  
        h.pt3dadd(-8510.041,5262.981,-1306.6105,self.nodeD, sec = self.node[0])

        h.pt3dadd(-8510.041,5262.981,-1306.6105,self.paraD1, sec = self.MYSA[0]) 
        h.pt3dadd(-8508.088,5261.5652,-1308.3941,self.paraD1, sec = self.MYSA[0])

        h.pt3dadd(-8508.088,5261.5652,-1308.3941,self.paraD2, sec = self.FLUT[0])
        h.pt3dadd(-8511.3668,5252.7289,-1311.7364,self.paraD2, sec = self.FLUT[0])

        h.pt3dadd(-8511.3668,5252.7289,-1311.7364,self.axonD, sec = self.STIN[0])
        h.pt3dadd(-8564.2095,5241.4464,-1332.816,self.axonD, sec = self.STIN[0])

        h.pt3dadd(-8564.2095,5241.4464,-1332.816,self.axonD, sec = self.STIN[1]) 
        h.pt3dadd(-8619.8514,5257.5034,-1336.0015,self.axonD, sec = self.STIN[1])

        h.pt3dadd(-8619.8514,5257.5034,-1336.0015,self.axonD, sec = self.STIN[2])
        h.pt3dadd(-8666.614,5262.5054,-1302.0569,self.axonD, sec = self.STIN[2])

        h.pt3dadd(-8666.614,5262.5054,-1302.0569,self.paraD2, sec = self.FLUT[1]) 
        h.pt3dadd(-8674.6766,5263.3678,-1296.2044,self.paraD2, sec = self.FLUT[1])

        h.pt3dadd(-8674.6766,5263.3678,-1296.2044,self.paraD1, sec = self.MYSA[1])
        h.pt3dadd(-8677.0953,5263.6265,-1294.4486,self.paraD1, sec = self.MYSA[1])

        h.pt3dadd(-8677.0953,5263.6265,-1294.4486,self.nodeD, sec = self.node[1])
        h.pt3dadd(-8677.9016,5263.7127,-1293.8634,self.nodeD, sec = self.node[1])

        h.pt3dadd(-8677.9016,5263.7127,-1293.8634,self.paraD1, sec = self.MYSA[2])
        h.pt3dadd(-8680.3203,5263.9714,-1292.1076,self.paraD1, sec = self.MYSA[2])

        h.pt3dadd(-8680.3203,5263.9714,-1292.1076,self.paraD2, sec = self.FLUT[2])
        h.pt3dadd(-8688.3828,5264.8338,-1286.2551,self.paraD2, sec = self.FLUT[2])

        h.pt3dadd(-8688.3828,5264.8338,-1286.2551,self.axonD, sec = self.STIN[3])
        h.pt3dadd(-8664.6489,5259.123,-1233.6425,self.axonD, sec = self.STIN[3])

        h.pt3dadd(-8664.6489,5259.123,-1233.6425,self.axonD, sec = self.STIN[4])
        h.pt3dadd(-8640.9149,5253.4121,-1181.0299,self.axonD, sec = self.STIN[4])

        h.pt3dadd(-8640.9149,5253.4121,-1181.0299,self.axonD, sec = self.STIN[5])
        h.pt3dadd(-8654.4099,5252.3139,-1124.6324,self.axonD, sec = self.STIN[5])

        h.pt3dadd(-8654.4099,5252.3139,-1124.6324,self.paraD2, sec = self.FLUT[3])
        h.pt3dadd(-8659.8731,5253.0038,-1116.2851,self.paraD2, sec = self.FLUT[3])

        h.pt3dadd(-8659.8731,5253.0038,-1116.2851,self.paraD1, sec = self.MYSA[3])
        h.pt3dadd(-8661.5121,5253.2107,-1113.7809,self.paraD1, sec = self.MYSA[3])

        h.pt3dadd(-8661.5121,5253.2107,-1113.7809,self.nodeD, sec = self.node[2])
        h.pt3dadd(-8662.0584,5253.2797,-1112.9462,self.nodeD, sec = self.node[2])

        h.pt3dadd(-8662.0584,5253.2797,-1112.9462,self.paraD1, sec = self.MYSA[4])
        h.pt3dadd(-8663.6973,5253.4866,-1110.442,self.paraD1, sec = self.MYSA[4])

        h.pt3dadd(-8663.6973,5253.4866,-1110.442,self.paraD2, sec = self.FLUT[4])
        h.pt3dadd(-8669.1605,5254.1765,-1102.0947,self.paraD2, sec = self.FLUT[4])

        h.pt3dadd(-8669.1605,5254.1765,-1102.0947,self.axonD, sec = self.STIN[6])
        h.pt3dadd(-8689.4073,5227.446,-1054.7708,self.axonD, sec = self.STIN[6])

        h.pt3dadd(-8689.4073,5227.446,-1054.7708,self.axonD, sec = self.STIN[7])
        h.pt3dadd(-8715.6231,5209.5333,-1006.2335,self.axonD, sec = self.STIN[7])

        h.pt3dadd(-8715.6231,5209.5333,-1006.2335,self.axonD, sec = self.STIN[8])
        h.pt3dadd(-8735.5594,5168.6021,-970.3005,self.axonD, sec = self.STIN[8])

        h.pt3dadd(-8735.5594,5168.6021,-970.3005,self.paraD2, sec = self.FLUT[5])
        h.pt3dadd(-8739.6336,5162.942,-963.1337,self.paraD2, sec = self.FLUT[5])

        h.pt3dadd(-8739.6336,5162.942,-963.1337,self.paraD1, sec = self.MYSA[5])
        h.pt3dadd(-8740.8559,5161.2439,-960.9836,self.paraD1, sec = self.MYSA[5])

        h.pt3dadd(-8740.8559,5161.2439,-960.9836,self.nodeD, sec = self.node[3])
        h.pt3dadd(-8741.2633,5160.6779,-960.2669,self.nodeD, sec = self.node[3])

        h.pt3dadd(-8741.2633,5160.6779,-960.2669,self.paraD1, sec = self.MYSA[6])
        h.pt3dadd(-8742.4856,5158.9799,-958.1169,self.paraD1, sec = self.MYSA[6])

        h.pt3dadd(-8742.4856,5158.9799,-958.1169,self.paraD2, sec = self.FLUT[6])
        h.pt3dadd(-8746.5599,5153.3197,-950.9501,self.paraD2, sec = self.FLUT[6])

        h.pt3dadd(-8746.5599,5153.3197,-950.9501,self.axonD, sec = self.STIN[9])
        h.pt3dadd(-8756.0258,5104.0933,-921.7755,self.axonD, sec = self.STIN[9])

        h.pt3dadd(-8756.0258,5104.0933,-921.7755,self.axonD, sec = self.STIN[10])
        h.pt3dadd(-8746.4489,5047.3709,-914.3698,self.axonD, sec = self.STIN[10])

        h.pt3dadd(-8746.4489,5047.3709,-914.3698,self.axonD, sec = self.STIN[11])
        h.pt3dadd(-8737.6291,4991.0199,-924.8948,self.axonD, sec = self.STIN[11])

        h.pt3dadd(-8737.6291,4991.0199,-924.8948,self.paraD2, sec = self.FLUT[7])
        h.pt3dadd(-8736.1085,4981.3042,-926.7095,self.paraD2, sec = self.FLUT[7])

        h.pt3dadd(-8736.1085,4981.3042,-926.7095,self.paraD1, sec = self.MYSA[7])
        h.pt3dadd(-8735.6523,4978.3895,-927.2539,self.paraD1, sec = self.MYSA[7])

        h.pt3dadd(-8735.6523,4978.3895,-927.2539,self.nodeD, sec = self.node[4])
        h.pt3dadd(-8735.5002,4977.4179,-927.4354,self.nodeD, sec = self.node[4])

        h.pt3dadd(-8735.5002,4977.4179,-927.4354,self.paraD1, sec = self.MYSA[8])
        h.pt3dadd(-8735.044,4974.5032,-927.9798,self.paraD1, sec = self.MYSA[8])

        h.pt3dadd(-8735.044,4974.5032,-927.9798,self.paraD2, sec = self.FLUT[8])
        h.pt3dadd(-8733.5234,4964.7875,-929.7944,self.paraD2, sec = self.FLUT[8])

        h.pt3dadd(-8733.5234,4964.7875,-929.7944,self.axonD, sec = self.STIN[12])
        h.pt3dadd(-8695.0828,4931.907,-958.1703,self.axonD, sec = self.STIN[12])

        h.pt3dadd(-8695.0828,4931.907,-958.1703,self.axonD, sec = self.STIN[13])
        h.pt3dadd(-8656.2785,4919.5846,-999.4787,self.axonD, sec = self.STIN[13])
        
        h.pt3dadd(-8656.2785,4919.5846,-999.4787,self.axonD, sec = self.STIN[14])
        h.pt3dadd(-8631.8627,4928.7037,-1051.2929,self.axonD, sec = self.STIN[14])
        
        h.pt3dadd(-8631.8627,4928.7037,-1051.2929,self.paraD2, sec = self.FLUT[9])
        h.pt3dadd(-8635.1668,4934.0557,-1059.0672,self.paraD2, sec = self.FLUT[9])
        
        h.pt3dadd(-8635.1668,4934.0557,-1059.0672,self.paraD1, sec = self.MYSA[9])
        h.pt3dadd(-8635.1806,4936.5603,-1060.7185,self.paraD1, sec = self.MYSA[9])
        
        h.pt3dadd(-8635.1806,4936.5603,-1060.7185,self.nodeD, sec = self.node[5])
        h.pt3dadd(-8635.1853,4937.3951,-1061.269,self.nodeD, sec = self.node[5])
        
        h.pt3dadd(-8635.1853,4937.3951,-1061.269,self.paraD1, sec = self.MYSA[10])
        h.pt3dadd(-8635.1991,4939.8997,-1062.9203,self.paraD1, sec = self.MYSA[10])
        
        h.pt3dadd(-8635.1991,4939.8997,-1062.9203,self.paraD2, sec = self.FLUT[10])
        h.pt3dadd(-8635.2454,4948.2484,-1068.4247,self.paraD2, sec = self.FLUT[10])
        
        h.pt3dadd(-8635.2454,4948.2484,-1068.4247,self.axonD, sec = self.STIN[15])
        h.pt3dadd(-8594.6456,4948.2792,-1109.8451,self.axonD, sec = self.STIN[15])
       
        h.pt3dadd(-8594.6456,4948.2792,-1109.8451,self.axonD, sec = self.STIN[16])
        h.pt3dadd(-8616.1496,4924.6866,-1158.27,self.axonD, sec = self.STIN[16])
        
        h.pt3dadd(-8616.1496,4924.6866,-1158.27,self.axonD, sec = self.STIN[17])
        h.pt3dadd(-8622.0889,4900.0709,-1210.4504,self.axonD, sec = self.STIN[17])

        h.pt3dadd(-8622.0889,4900.0709,-1210.4504,self.paraD2, sec = self.FLUT[11])
        h.pt3dadd(-8614.7535,4893.327,-1209.6066,self.paraD2, sec = self.FLUT[11])

        h.pt3dadd(-8614.7535,4893.327,-1209.6066,self.paraD1, sec = self.MYSA[11])
        h.pt3dadd(-8612.5529,4891.3039,-1209.3534,self.paraD1, sec = self.MYSA[11])

        h.pt3dadd(-8612.5529,4891.3039,-1209.3534,self.nodeD, sec = self.node[6])
        h.pt3dadd(-8611.8974,4890.5864,-1209.5891,self.nodeD, sec = self.node[6])

        h.pt3dadd(-8611.8974,4890.5864,-1209.5891,self.paraD1, sec = self.MYSA[12])
        h.pt3dadd(-8609.0417,4889.6727,-1209.4864,self.paraD1, sec = self.MYSA[12])

        h.pt3dadd(-8609.0417,4889.6727,-1209.4864,self.paraD2, sec = self.FLUT[12])
        h.pt3dadd(-8599.5229,4886.6272,-1209.1441,self.paraD2, sec = self.FLUT[12])

        h.pt3dadd(-8599.5229,4886.6272,-1209.1441,self.axonD, sec = self.STIN[18])
        h.pt3dadd(-8548.6859,4858.7188,-1208.299,self.axonD, sec = self.STIN[18])

        h.pt3dadd(-8548.6859,4858.7188,-1208.299,self.axonD, sec = self.STIN[19])
        h.pt3dadd(-8496.2132,4834.803,-1214.5164,self.axonD, sec = self.STIN[19])

        h.pt3dadd(-8496.2132,4834.803,-1214.5164,self.axonD, sec = self.STIN[20])
        h.pt3dadd(-8441.6614,4815.1775,-1216.2317,self.axonD, sec = self.STIN[20])

        h.pt3dadd(-8441.6614,4815.1775,-1216.2317,self.paraD2, sec = self.FLUT[13])
        h.pt3dadd(-8434.2338,4811.2573,-1210.8039,self.paraD2, sec = self.FLUT[13])

        h.pt3dadd(-8434.2338,4811.2573,-1210.8039,self.paraD1, sec = self.MYSA[13])
        h.pt3dadd(-8432.0055,4810.0813,-1209.1755,self.paraD1, sec = self.MYSA[13])

        h.pt3dadd(-8432.0055,4810.0813,-1209.1755,self.nodeD, sec = self.node[7])
        h.pt3dadd(-8431.2627,4809.6892,-1208.6327,self.nodeD, sec = self.node[7])

        h.pt3dadd(-8431.2627,4809.6892,-1208.6327,self.paraD1, sec = self.MYSA[14])
        h.pt3dadd(-8429.0344,4808.5132,-1207.0044,self.paraD1, sec = self.MYSA[14])

        h.pt3dadd(-8429.0344,4808.5132,-1207.0044,self.paraD2, sec = self.FLUT[14])
        h.pt3dadd(-8421.6067,4804.5929,-1201.5765,self.paraD2, sec = self.FLUT[14])

        h.pt3dadd(-8421.6067,4804.5929,-1201.5765,self.axonD, sec = self.STIN[21])
        h.pt3dadd(-8370.2149,4779.4358,-1192.09,self.axonD, sec = self.STIN[21])

        h.pt3dadd(-8370.2149,4779.4358,-1192.09,self.axonD, sec = self.STIN[22])
        h.pt3dadd(-8318.8562,4755.2471,-1180.2076,self.axonD, sec = self.STIN[22])

        h.pt3dadd(-8318.8562,4755.2471,-1180.2076,self.axonD, sec = self.STIN[23])
        h.pt3dadd(-8275.0906,4717.9033,-1187.5567,self.axonD, sec = self.STIN[23])

        h.pt3dadd(-8275.0906,4717.9033,-1187.5567,self.paraD2, sec = self.FLUT[15])
        h.pt3dadd(-8265.7935,4715.1591,-1185.1004,self.paraD2, sec = self.FLUT[15])

        h.pt3dadd(-8265.7935,4715.1591,-1185.1004,self.paraD1, sec = self.MYSA[15])
        h.pt3dadd(-8263.0043,4714.3358,-1184.3636,self.paraD1, sec = self.MYSA[15])

        h.pt3dadd(-8263.0043,4714.3358,-1184.3636,self.nodeD, sec = self.node[8])
        h.pt3dadd(-8262.0746,4714.0614,-1184.1179,self.nodeD, sec = self.node[8])

        h.pt3dadd(-8262.0746,4714.0614,-1184.1179,self.paraD1, sec = self.MYSA[16])
        h.pt3dadd(-8259.2855,4713.2381,-1183.3811,self.paraD1, sec = self.MYSA[16])

        h.pt3dadd(-8259.2855,4713.2381,-1183.3811,self.paraD2, sec = self.FLUT[16])
        h.pt3dadd(-8249.501,4711.1804,-1183.2106,self.paraD2, sec = self.FLUT[16])

        h.pt3dadd(-8249.501,4711.1804,-1183.2106,self.axonD, sec = self.STIN[24])
        h.pt3dadd(-8196.5028,4689.7858,-1173.3385,self.axonD, sec = self.STIN[24])

        h.pt3dadd(-8196.5028,4689.7858,-1173.3385,self.axonD, sec = self.STIN[25])
        h.pt3dadd(-8139.6545,4684.6762,-1163.0347,self.axonD, sec = self.STIN[25])

        h.pt3dadd(-8139.6545,4684.6762,-1163.0347,self.axonD, sec = self.STIN[26])
        h.pt3dadd(-8083.1797,4671.5268,-1164.3368,self.axonD, sec = self.STIN[26])

        h.pt3dadd(-8083.1797,4671.5268,-1164.3368,self.paraD2, sec = self.FLUT[17])
        h.pt3dadd(-8073.4427,4669.2596,-1164.5613,self.paraD2, sec = self.FLUT[17])

        h.pt3dadd(-8073.4427,4669.2596,-1164.5613,self.paraD1, sec = self.MYSA[17])
        h.pt3dadd(-8070.5216,4668.5795,-1164.6286,self.paraD1, sec = self.MYSA[17])

        h.pt3dadd(-8070.5216,4668.5795,-1164.6286,self.nodeD, sec = self.node[9])
        h.pt3dadd(-8069.5479,4668.3528,-1164.6511,self.nodeD, sec = self.node[9])

        h.pt3dadd(-8069.5479,4668.3528,-1164.6511,self.paraD1, sec = self.MYSA[18])
        h.pt3dadd(-8066.7504,4667.3023,-1164.3862,self.paraD1, sec = self.MYSA[18])

        h.pt3dadd(-8066.7504,4667.3023,-1164.3862,self.paraD2, sec = self.FLUT[18])
        h.pt3dadd(-8057.4251,4663.8007,-1163.5033,self.paraD2, sec = self.FLUT[18])

        h.pt3dadd(-8057.4251,4663.8007,-1163.5033,self.axonD, sec = self.STIN[27])
        h.pt3dadd(-8014.4343,4653.8019,-1201.1305,self.axonD, sec = self.STIN[27])

        h.pt3dadd(-8014.4343,4653.8019,-1201.1305,self.axonD, sec = self.STIN[28])
        h.pt3dadd(-7959.7769,4642.2735,-1216.7403,self.axonD, sec = self.STIN[28])

        h.pt3dadd(-7959.7769,4642.2735,-1216.7403,self.axonD, sec = self.STIN[29])
        h.pt3dadd(-7904.8252,4626.7471,-1206.5794,self.axonD, sec = self.STIN[29])

        h.pt3dadd(-7904.8252,4626.7471,-1206.5794,self.paraD2, sec = self.FLUT[19])
        h.pt3dadd(-7895.3508,4624.0702,-1204.8275,self.paraD2, sec = self.FLUT[19])

        h.pt3dadd(-7895.3508,4624.0702,-1204.8275,self.paraD1, sec = self.MYSA[19])
        h.pt3dadd(-7892.5418,4623.053,-1204.5537,self.paraD1, sec = self.MYSA[19])

        h.pt3dadd(-7892.5418,4623.053,-1204.5537,self.nodeD, sec = self.node[10])
        h.pt3dadd(-7891.6055,4622.7139,-1204.4624,self.nodeD, sec = self.node[10])

        h.pt3dadd(-7891.6055,4622.7139,-1204.4624,self.paraD1, sec = self.MYSA[20])
        h.pt3dadd(-7888.7965,4621.6967,-1204.1886,self.paraD1, sec = self.MYSA[20])

        h.pt3dadd(-7888.7965,4621.6967,-1204.1886,self.paraD2, sec = self.FLUT[20])
        h.pt3dadd(-7880.3929,4616.4506,-1202.8258,self.paraD2, sec = self.FLUT[20])

        h.pt3dadd(-7880.3929,4616.4506,-1202.8258,self.axonD, sec = self.STIN[30])
        h.pt3dadd(-7823.8781,4614.0013,-1190.0163,self.axonD, sec = self.STIN[30])

        h.pt3dadd(-7823.8781,4614.0013,-1190.0163,self.axonD, sec = self.STIN[31])
        h.pt3dadd(-7769.83,4603.5238,-1171.7674,self.axonD, sec = self.STIN[31])

        h.pt3dadd(-7769.83,4603.5238,-1171.7674,self.axonD, sec = self.STIN[32])
        h.pt3dadd(-7719.2153,4576.3835,-1163.6706,self.axonD, sec = self.STIN[32])

        h.pt3dadd(-7719.2153,4576.3835,-1163.6706,self.paraD2, sec = self.FLUT[21])
        h.pt3dadd(-7710.4886,4571.7042,-1162.2746,self.paraD2, sec = self.FLUT[21])

        h.pt3dadd(-7710.4886,4571.7042,-1162.2746,self.paraD1, sec = self.MYSA[21])
        h.pt3dadd(-7707.8706,4570.3003,-1161.8558,self.paraD1, sec = self.MYSA[21])

        h.pt3dadd(-7707.8706,4570.3003,-1161.8558,self.nodeD, sec = self.node[11])
        h.pt3dadd(-7706.998,4569.8324,-1161.7162,self.nodeD, sec = self.node[11])

        h.pt3dadd(-7706.998,4569.8324,-1161.7162,self.paraD1, sec = self.MYSA[22])
        h.pt3dadd(-7704.4927,4571.1495,-1160.7218,self.paraD1, sec = self.MYSA[22])

        h.pt3dadd(-7704.4927,4571.1495,-1160.7218,self.paraD2, sec = self.FLUT[22])
        h.pt3dadd(-7696.1419,4575.54,-1157.407,self.paraD2, sec = self.FLUT[22])

        h.pt3dadd(-7696.1419,4575.54,-1157.407,self.axonD, sec = self.STIN[33])
        h.pt3dadd(-7647.2578,4560.6046,-1129.9975,self.axonD, sec = self.STIN[33])

        h.pt3dadd(-7647.2578,4560.6046,-1129.9975,self.axonD, sec = self.STIN[34])
        h.pt3dadd(-7599.2942,4534.0076,-1111.1273,self.axonD, sec = self.STIN[34])

        h.pt3dadd(-7599.2942,4534.0076,-1111.1273,self.axonD, sec = self.STIN[35])
        h.pt3dadd(-7554.2687,4513.5498,-1080.8256,self.axonD, sec = self.STIN[35])

        h.pt3dadd(-7554.2687,4513.5498,-1080.8256,self.paraD2, sec = self.FLUT[23])
        h.pt3dadd(-7548.6069,4507.6695,-1075.0494,self.paraD2, sec = self.FLUT[23])

        h.pt3dadd(-7548.6069,4507.6695,-1075.0494,self.paraD1, sec = self.MYSA[23])
        h.pt3dadd(-7546.9083,4505.9054,-1073.3165,self.paraD1, sec = self.MYSA[23])

        h.pt3dadd(-7546.9083,4505.9054,-1073.3165,self.nodeD, sec = self.node[12])
        h.pt3dadd(-7546.3421,4505.3173,-1072.7389,self.nodeD, sec = self.node[12])

        h.pt3dadd(-7546.3421,4505.3173,-1072.7389,self.paraD1, sec = self.MYSA[24])
        h.pt3dadd(-7543.63,4505.9373,-1071.6162,self.paraD1, sec = self.MYSA[24])

        h.pt3dadd(-7543.63,4505.9373,-1071.6162,self.paraD2, sec = self.FLUT[24])
        h.pt3dadd(-7533.8449,4505.1767,-1069.6998,self.paraD2, sec = self.FLUT[24])

        h.pt3dadd(-7533.8449,4505.1767,-1069.6998,self.axonD, sec = self.STIN[36])
        h.pt3dadd(-7477.0443,4494.8706,-1064.0891,self.axonD, sec = self.STIN[36])

        h.pt3dadd(-7477.0443,4494.8706,-1064.0891,self.axonD, sec = self.STIN[37])
        h.pt3dadd(-7421.0341,4492.0227,-1078.8791,self.axonD, sec = self.STIN[37])

        h.pt3dadd(-7421.0341,4492.0227,-1078.8791,self.axonD, sec = self.STIN[38])
        h.pt3dadd(-7380.8116,4528.3714,-1099.4928,self.axonD, sec = self.STIN[38])

        h.pt3dadd(-7380.8116,4528.3714,-1099.4928,self.paraD2, sec = self.FLUT[25])
        h.pt3dadd(-7371.9149,4527.8995,-1094.9513,self.paraD2, sec = self.FLUT[25])

        h.pt3dadd(-7371.9149,4527.8995,-1094.9513,self.paraD1, sec = self.MYSA[25])
        h.pt3dadd(-7369.2459,4527.758,-1093.5888,self.paraD1, sec = self.MYSA[25])

        h.pt3dadd(-7369.2459,4527.758,-1093.5888,self.nodeD, sec = self.node[13])
        h.pt3dadd(-7368.3562,4527.7108,-1093.1347,self.nodeD, sec = self.node[13])

        h.pt3dadd(-7368.3562,4527.7108,-1093.1347,self.paraD1, sec = self.MYSA[26])
        h.pt3dadd(-7365.6872,4527.5692,-1091.7722,self.paraD1, sec = self.MYSA[26])

        h.pt3dadd(-7365.6872,4527.5692,-1091.7722,self.paraD2, sec = self.FLUT[26])
        h.pt3dadd(-7355.7274,4526.872,-1091.2093,self.paraD2, sec = self.FLUT[26])

        h.pt3dadd(-7355.7274,4526.872,-1091.2093,self.axonD, sec = self.STIN[39])
        h.pt3dadd(-7314.9684,4520.8258,-1050.3909,self.axonD, sec = self.STIN[39])

        h.pt3dadd(-7314.9684,4520.8258,-1050.3909,self.axonD, sec = self.STIN[40])
        h.pt3dadd(-7292.2689,4493.7285,-1004.4075,self.axonD, sec = self.STIN[40])

        h.pt3dadd(-7292.2689,4493.7285,-1004.4075,self.axonD, sec = self.STIN[41])
        h.pt3dadd(-7257.9313,4484.1865,-958.6486,self.axonD, sec = self.STIN[41])

        h.pt3dadd(-7257.9313,4484.1865,-958.6486,self.paraD2, sec = self.FLUT[27])
        h.pt3dadd(-7251.3663,4478.6433,-953.5325,self.paraD2, sec = self.FLUT[27])

        h.pt3dadd(-7251.3663,4478.6433,-953.5325,self.paraD1, sec = self.MYSA[27])
        h.pt3dadd(-7249.3968,4476.9803,-951.9977,self.paraD1, sec = self.MYSA[27])

        h.pt3dadd(-7249.3968,4476.9803,-951.9977,self.nodeD, sec = self.node[14])
        h.pt3dadd(-7248.7403,4476.426,-951.4861,self.nodeD, sec = self.node[14])

        h.pt3dadd(-7248.7403,4476.426,-951.4861,self.paraD1, sec = self.MYSA[28])
        h.pt3dadd(-7247.8585,4475.0767,-948.9559,self.paraD1, sec = self.MYSA[28])

        h.pt3dadd(-7247.8585,4475.0767,-948.9559,self.paraD2, sec = self.FLUT[28])
        h.pt3dadd(-7240.7511,4473.1346,-942.1947,self.paraD2, sec = self.FLUT[28])

        h.pt3dadd(-7240.7511,4473.1346,-942.1947,self.axonD, sec = self.STIN[42])
        h.pt3dadd(-7215.6747,4459.4884,-891.7075,self.axonD, sec = self.STIN[42])

        h.pt3dadd(-7215.6747,4459.4884,-891.7075,self.axonD, sec = self.STIN[43])
        h.pt3dadd(-7197.7272,4472.545,-838.122,self.axonD, sec = self.STIN[43])

        h.pt3dadd(-7197.7272,4472.545,-838.122,self.axonD, sec = self.STIN[44])
        h.pt3dadd(-7180.3048,4491.4892,-786.1453,self.axonD, sec = self.STIN[44])

        h.pt3dadd(-7180.3048,4491.4892,-786.1453,self.paraD2, sec = self.FLUT[29])
        h.pt3dadd(-7172.0616,4494.2692,-781.2136,self.paraD2, sec = self.FLUT[29])

        h.pt3dadd(-7172.0616,4494.2692,-781.2136,self.paraD1, sec = self.MYSA[29])
        h.pt3dadd(-7169.5886,4495.1032,-779.7341,self.paraD1, sec = self.MYSA[29])

        h.pt3dadd(-7169.5886,4495.1032,-779.7341,self.nodeD, sec = self.node[15])
        h.pt3dadd(-7168.7643,4495.3812,-779.241,self.nodeD, sec = self.node[15])

        h.pt3dadd(-7168.7643,4495.3812,-779.241,self.paraD1, sec = self.MYSA[30])
        h.pt3dadd(-7166.2913,4496.2152,-777.7615,self.paraD1, sec = self.MYSA[30])

        h.pt3dadd(-7166.2913,4496.2152,-777.7615,self.paraD2, sec = self.FLUT[30])
        h.pt3dadd(-7158.0481,4498.9951,-772.8298,self.paraD2, sec = self.FLUT[30])

        h.pt3dadd(-7158.0481,4498.9951,-772.8298,self.axonD, sec = self.STIN[45])
        h.pt3dadd(-7122.9562,4522.1368,-732.867,self.axonD, sec = self.STIN[45])

        h.pt3dadd(-7122.9562,4522.1368,-732.867,self.axonD, sec = self.STIN[46])
        h.pt3dadd(-7103.4585,4552.8665,-687.706,self.axonD, sec = self.STIN[46])

        h.pt3dadd(-7103.4585,4552.8665,-687.706,self.axonD, sec = self.STIN[47])
        h.pt3dadd(-7058.9961,4555.8978,-650.5858,self.axonD, sec = self.STIN[47])

        h.pt3dadd(-7058.9961,4555.8978,-650.5858,self.paraD2, sec = self.FLUT[31])
        h.pt3dadd(-7051.3302,4556.4204,-644.1858,self.paraD2, sec = self.FLUT[31])

        h.pt3dadd(-7051.3302,4556.4204,-644.1858,self.paraD1, sec = self.MYSA[31])
        h.pt3dadd(-7049.0054,4556.9396,-642.362,self.paraD1, sec = self.MYSA[31])

        h.pt3dadd(-7049.0054,4556.9396,-642.362,self.nodeD, sec = self.node[16])
        h.pt3dadd(-7048.2305,4557.1127,-641.7541,self.nodeD, sec = self.node[16])

        h.pt3dadd(-7048.2305,4557.1127,-641.7541,self.paraD1, sec = self.MYSA[32])
        h.pt3dadd(-7049.7971,4556.424,-639.2901,self.paraD1, sec = self.MYSA[32])

        h.pt3dadd(-7049.7971,4556.424,-639.2901,self.paraD2, sec = self.FLUT[32])
        h.pt3dadd(-7055.0189,4554.1283,-631.0765,self.paraD2, sec = self.FLUT[32])

        h.pt3dadd(-7055.0189,4554.1283,-631.0765,self.axonD, sec = self.STIN[48])
        h.pt3dadd(-7028.9682,4552.8519,-579.2717,self.axonD, sec = self.STIN[48])

        h.pt3dadd(-7028.9682,4552.8519,-579.2717,self.axonD, sec = self.STIN[49])
        h.pt3dadd(-6996.2317,4567.7239,-533.7619,self.axonD, sec = self.STIN[49])

        h.pt3dadd(-6996.2317,4567.7239,-533.7619,self.axonD, sec = self.STIN[50])
        h.pt3dadd(-6980.168,4588.9325,-482.224,self.axonD, sec = self.STIN[50])

        h.pt3dadd(-6980.168,4588.9325,-482.224,self.paraD2, sec = self.FLUT[33])
        h.pt3dadd(-6975.542,4595.8564,-476.6869,self.paraD2, sec = self.FLUT[33])

        h.pt3dadd(-6975.542,4595.8564,-476.6869,self.paraD1, sec = self.MYSA[33])
        h.pt3dadd(-6976.0926,4595.4393,-473.7675,self.paraD1, sec = self.MYSA[33])

        h.pt3dadd(-6976.0926,4595.4393,-473.7675,self.nodeD, sec = self.node[17])
        h.pt3dadd(-6976.2761,4595.3002,-472.7944,self.nodeD, sec = self.node[17])

        h.pt3dadd(-6976.2761,4595.3002,-472.7944,self.paraD1, sec = self.MYSA[34])
        h.pt3dadd(-6976.6384,4596.592,-470.111,self.paraD1, sec = self.MYSA[34])

        h.pt3dadd(-6976.6384,4596.592,-470.111,self.paraD2, sec = self.FLUT[34])
        h.pt3dadd(-6977.846,4600.8976,-461.1666,self.paraD2, sec = self.FLUT[34])

        h.pt3dadd(-6977.846,4600.8976,-461.1666,self.axonD, sec = self.STIN[51])
        h.pt3dadd(-6954.8127,4607.4114,-408.3363,self.axonD, sec = self.STIN[51])

        h.pt3dadd(-6954.8127,4607.4114,-408.3363,self.axonD, sec = self.STIN[52])
        h.pt3dadd(-6928.6666,4605.0873,-356.6161,self.axonD, sec = self.STIN[52])

        h.pt3dadd(-6928.6666,4605.0873,-356.6161,self.axonD, sec = self.STIN[53])
        h.pt3dadd(-6926.5936,4609.8572,-298.8498,self.axonD, sec = self.STIN[53])

        h.pt3dadd(-6926.5936,4609.8572,-298.8498,self.paraD2, sec = self.FLUT[35])
        h.pt3dadd(-6926.2362,4610.6796,-288.89,self.paraD2, sec = self.FLUT[35])

        h.pt3dadd(-6926.2362,4610.6796,-288.89,self.paraD1, sec = self.MYSA[35])
        h.pt3dadd(-6926.129,4610.9264,-285.9021,self.paraD1, sec = self.MYSA[35])

        h.pt3dadd(-6926.129,4610.9264,-285.9021,self.nodeD, sec = self.node[18])
        h.pt3dadd(-6926.0933,4611.0086,-284.9062,self.nodeD, sec = self.node[18])

        h.pt3dadd(-6926.0933,4611.0086,-284.9062,self.paraD1, sec = self.MYSA[36])
        h.pt3dadd(-6925.986,4611.2553,-281.9182,self.paraD1, sec = self.MYSA[36])

        h.pt3dadd(-6925.986,4611.2553,-281.9182,self.paraD2, sec = self.FLUT[36])
        h.pt3dadd(-6925.6286,4612.0777,-271.9585,self.paraD2, sec = self.FLUT[36])

        h.pt3dadd(-6925.6286,4612.0777,-271.9585,self.axonD, sec = self.STIN[54])
        h.pt3dadd(-6923.4933,4600.9334,-215.0793,self.axonD, sec = self.STIN[54])

        h.pt3dadd(-6923.4933,4600.9334,-215.0793,self.axonD, sec = self.STIN[55])
        h.pt3dadd(-6907.9839,4594.5392,-159.5584,self.axonD, sec = self.STIN[55])

        h.pt3dadd(-6907.9839,4594.5392,-159.5584,self.axonD, sec = self.STIN[56])
        h.pt3dadd(-6865.3449,4609.5922,-123.2358,self.axonD, sec = self.STIN[56])

        h.pt3dadd(-6865.3449,4609.5922,-123.2358,self.paraD2, sec = self.FLUT[37])
        h.pt3dadd(-6859.59,4602.7279,-118.7902,self.paraD2, sec = self.FLUT[37])

        h.pt3dadd(-6859.59,4602.7279,-118.7902,self.paraD1, sec = self.MYSA[37])
        h.pt3dadd(-6856.6299,4602.4873,-119.2145,self.paraD1, sec = self.MYSA[37])

        h.pt3dadd(-6856.6299,4602.4873,-119.2145,self.nodeD, sec = self.node[19])
        h.pt3dadd(-6855.6432,4602.407,-119.3559,self.nodeD, sec = self.node[19])

        h.pt3dadd(-6855.6432,4602.407,-119.3559,self.paraD1, sec = self.MYSA[38])
        h.pt3dadd(-6852.6831,4602.1663,-119.7801,self.paraD1, sec = self.MYSA[38])

        h.pt3dadd(-6852.6831,4602.1663,-119.7801,self.paraD2, sec = self.FLUT[38])
        h.pt3dadd(-6846.5381,4595.9263,-114.953,self.paraD2, sec = self.FLUT[38])

        h.pt3dadd(-6846.5381,4595.9263,-114.953,self.axonD, sec = self.STIN[57])
        h.pt3dadd(-6830.3341,4588.0068,-59.8285,self.axonD, sec = self.STIN[57])

        h.pt3dadd(-6830.3341,4588.0068,-59.8285,self.axonD, sec = self.STIN[58])
        h.pt3dadd(-6798.0419,4594.0458,-12.0295,self.axonD, sec = self.STIN[58])

        h.pt3dadd(-6798.0419,4594.0458,-12.0295,self.axonD, sec = self.STIN[59])
        h.pt3dadd(-6794.9846,4586.6646,45.4177,self.axonD, sec = self.STIN[59])

        h.pt3dadd(-6794.9846,4586.6646,45.4177,self.paraD2, sec = self.FLUT[39])
        h.pt3dadd(-6795.0024,4584.603,55.2028,self.paraD2, sec = self.FLUT[39])

        h.pt3dadd(-6795.0024,4584.603,55.2028,self.paraD1, sec = self.MYSA[39])
        h.pt3dadd(-6795.0078,4583.9845,58.1384,self.paraD1, sec = self.MYSA[39])

        h.pt3dadd(-6795.0078,4583.9845,58.1384,self.nodeD, sec = self.node[20])
        h.pt3dadd(-6795.0096,4583.7784,59.1169,self.nodeD, sec = self.node[20])

        h.pt3dadd(-6795.0096,4583.7784,59.1169,self.paraD1, sec = self.MYSA[40])
        h.pt3dadd(-6795.0149,4583.1599,62.0524,self.paraD1, sec = self.MYSA[40])

        h.pt3dadd(-6795.0149,4583.1599,62.0524,self.paraD2, sec = self.FLUT[40])
        h.pt3dadd(-6795.0328,4581.0983,71.8376,self.paraD2, sec = self.FLUT[40])

        h.pt3dadd(-6795.0328,4581.0983,71.8376,self.axonD, sec = self.STIN[60])
        h.pt3dadd(-6798.3696,4558.6999,125.234,self.axonD, sec = self.STIN[60])

        h.pt3dadd(-6798.3696,4558.6999,125.234,self.axonD, sec = self.STIN[61])
        h.pt3dadd(-6775.9139,4559.13,178.7089,self.axonD, sec = self.STIN[61])

        h.pt3dadd(-6775.9139,4559.13,178.7089,self.axonD, sec = self.STIN[62])
        h.pt3dadd(-6773.9046,4585.0023,230.5797,self.axonD, sec = self.STIN[62])

        h.pt3dadd(-6773.9046,4585.0023,230.5797,self.paraD2, sec = self.FLUT[41])
        h.pt3dadd(-6773.5582,4589.4631,239.523,self.paraD2, sec = self.FLUT[41])

        h.pt3dadd(-6773.5582,4589.4631,239.523,self.paraD1, sec = self.MYSA[41])
        h.pt3dadd(-6773.4543,4590.8013,242.2059,self.paraD1, sec = self.MYSA[41])

        h.pt3dadd(-6773.4543,4590.8013,242.2059,self.nodeD, sec = self.node[21])
        h.pt3dadd(-6773.4196,4591.2474,243.1003,self.nodeD, sec = self.node[21])

        h.pt3dadd(-6773.4196,4591.2474,243.1003,self.paraD1, sec = self.MYSA[42])
        h.pt3dadd(-6773.3157,4592.5856,245.7832,self.paraD1, sec = self.MYSA[42])

        h.pt3dadd(-6773.3157,4592.5856,245.7832,self.paraD2, sec = self.FLUT[42])
        h.pt3dadd(-6776.0921,4600.3697,251.4134,self.paraD2, sec = self.FLUT[42])

        h.pt3dadd(-6776.0921,4600.3697,251.4134,self.axonD, sec = self.STIN[63])
        h.pt3dadd(-6783.9602,4652.3206,275.9734,self.axonD, sec = self.STIN[63])

        h.pt3dadd(-6783.9602,4652.3206,275.9734,self.axonD, sec = self.STIN[64])
        h.pt3dadd(-6785.0654,4696.7541,313.2353,self.axonD, sec = self.STIN[64])

        h.pt3dadd(-6785.0654,4696.7541,313.2353,self.axonD, sec = self.STIN[65])
        h.pt3dadd(-6793.348,4747.11,340.7976,self.axonD, sec = self.STIN[65])

        h.pt3dadd(-6793.348,4747.11,340.7976,self.paraD2, sec = self.FLUT[43])
        h.pt3dadd(-6794.7761,4755.7921,345.5497,self.paraD2, sec = self.FLUT[43])

        h.pt3dadd(-6794.7761,4755.7921,345.5497,self.paraD1, sec = self.MYSA[43])
        h.pt3dadd(-6794.8653,4758.6706,346.3899,self.paraD1, sec = self.MYSA[43])
        
        h.pt3dadd(-6794.8653,4758.6706,346.3899,self.nodeD, sec = self.node[22])
        h.pt3dadd(-6794.8951,4759.6301,346.6699,self.nodeD, sec = self.node[22])
        
        h.pt3dadd(-6794.8951,4759.6301,346.6699,self.paraD1, sec = self.MYSA[44])
        h.pt3dadd(-6794.9844,4762.5087,347.5102,self.paraD1, sec = self.MYSA[44])
        
        h.pt3dadd(-6794.9844,4762.5087,347.5102,self.paraD2, sec = self.FLUT[44])
        h.pt3dadd(-6795.282,4772.1039,350.3109,self.paraD2, sec = self.FLUT[44])
        
        h.pt3dadd(-6795.282,4772.1039,350.3109,self.axonD, sec = self.STIN[66])
        h.pt3dadd(-6805.1904,4827.881,362.7502,self.axonD, sec = self.STIN[66])
       
        h.pt3dadd(-6805.1904,4827.881,362.7502,self.axonD, sec = self.STIN[67])
        h.pt3dadd(-6808.9733,4882.5417,381.7737,self.axonD, sec = self.STIN[67])
        
        h.pt3dadd(-6808.9733,4882.5417,381.7737,self.axonD, sec = self.STIN[68])
        h.pt3dadd(-6818.4085,4938.9449,391.4518,self.axonD, sec = self.STIN[68])
       
        h.pt3dadd(-6818.4085,4938.9449,391.4518,self.paraD2, sec = self.FLUT[45])
        h.pt3dadd(-6819.2908,4948.4522,394.4236,self.paraD2, sec = self.FLUT[45])
        
        h.pt3dadd(-6819.2908,4948.4522,394.4236,self.paraD1, sec = self.MYSA[45])
        h.pt3dadd(-6819.5555,4951.3045,395.3151,self.paraD1, sec = self.MYSA[45])
        
        h.pt3dadd(-6819.5555,4951.3045,395.3151,self.nodeD, sec = self.node[23])
        h.pt3dadd(-6819.6437,4952.2552,395.6122,self.nodeD, sec = self.node[23])

        h.pt3dadd(-6819.6437,4952.2552,395.6122,self.paraD1, sec = self.MYSA[46])
        h.pt3dadd(-6819.9084,4955.1074,396.5038,self.paraD1, sec = self.MYSA[46])

        h.pt3dadd(-6819.9084,4955.1074,396.5038,self.paraD2, sec = self.FLUT[46])
        h.pt3dadd(-6820.7907,4964.6148,399.4755,self.paraD2, sec = self.FLUT[46])

        h.pt3dadd(-6820.7907,4964.6148,399.4755,self.axonD, sec = self.STIN[69])
        h.pt3dadd(-6831.2312,5020.3109,411.8426,self.axonD, sec = self.STIN[69])

        h.pt3dadd(-6831.2312,5020.3109,411.8426,self.axonD, sec = self.STIN[70])
        h.pt3dadd(-6845.0752,5075.0888,424.9468,self.axonD, sec = self.STIN[70])

        h.pt3dadd(-6845.0752,5075.0888,424.9468,self.axonD, sec = self.STIN[71])
        h.pt3dadd(-6851.7066,5131.0344,438.7349,self.axonD, sec = self.STIN[71])

        h.pt3dadd(-6851.7066,5131.0344,438.7349,self.paraD2, sec = self.FLUT[47])
        h.pt3dadd(-6852.8499,5140.6802,441.1121,self.paraD2, sec = self.FLUT[47])

        h.pt3dadd(-6852.8499,5140.6802,441.1121,self.paraD1, sec = self.MYSA[47])
        h.pt3dadd(-6853.1929,5143.574,441.8253,self.paraD1, sec = self.MYSA[47])

        h.pt3dadd(-6853.1929,5143.574,441.8253,self.nodeD, sec = self.node[24])
        h.pt3dadd(-6853.3073,5144.5385,442.063,self.nodeD, sec = self.node[24])

        h.pt3dadd(-6853.3073,5144.5385,442.063,self.paraD1, sec = self.MYSA[48])
        h.pt3dadd(-6853.6503,5147.4323,442.7762,self.paraD1, sec = self.MYSA[48])

        h.pt3dadd(-6853.6503,5147.4323,442.7762,self.paraD2, sec = self.FLUT[48])
        h.pt3dadd(-6855.9157,5156.6918,445.7975,self.paraD2, sec = self.FLUT[48])

        h.pt3dadd(-6855.9157,5156.6918,445.7975,self.axonD, sec = self.STIN[72])
        h.pt3dadd(-6860.2423,5211.7482,463.5198,self.axonD, sec = self.STIN[72])

        h.pt3dadd(-6860.2423,5211.7482,463.5198,self.axonD, sec = self.STIN[73])
        h.pt3dadd(-6840.4626,5262.866,482.4865,self.axonD, sec = self.STIN[73])

        h.pt3dadd(-6840.4626,5262.866,482.4865,self.axonD, sec = self.STIN[74])
        h.pt3dadd(-6842.5123,5319.0834,496.6075,self.axonD, sec = self.STIN[74])

        h.pt3dadd(-6842.5123,5319.0834,496.6075,self.paraD2, sec = self.FLUT[49])
        h.pt3dadd(-6842.8657,5328.7761,499.0421,self.paraD2, sec = self.FLUT[49])

        h.pt3dadd(-6842.8657,5328.7761,499.0421,self.paraD1, sec = self.MYSA[49])
        h.pt3dadd(-6842.9717,5331.6838,499.7725,self.paraD1, sec = self.MYSA[49])

        h.pt3dadd(-6842.9717,5331.6838,499.7725,self.nodeD, sec = self.node[25])
        h.pt3dadd(-6843.0071,5332.6531,500.016,self.nodeD, sec = self.node[25])

        h.pt3dadd(-6843.0071,5332.6531,500.016,self.paraD1, sec = self.MYSA[50])
        h.pt3dadd(-6843.1131,5335.5609,500.7464,self.paraD1, sec = self.MYSA[50])

        h.pt3dadd(-6843.1131,5335.5609,500.7464,self.paraD2, sec = self.FLUT[50])
        h.pt3dadd(-6843.4665,5345.2536,503.181,self.paraD2, sec = self.FLUT[50])

        h.pt3dadd(-6843.4665,5345.2536,503.181,self.axonD, sec = self.STIN[75])
        h.pt3dadd(-6845.378,5403.0083,508.1552,self.axonD, sec = self.STIN[75])

        h.pt3dadd(-6845.378,5403.0083,508.1552,self.axonD, sec = self.STIN[76])
        h.pt3dadd(-6816.054,5452.918,504.5336,self.axonD, sec = self.STIN[76])

        h.pt3dadd(-6816.054,5452.918,504.5336,self.axonD, sec = self.STIN[77])
        h.pt3dadd(-6771.0747,5489.2932,508.7425,self.axonD, sec = self.STIN[77])

        h.pt3dadd(-6771.0747,5489.2932,508.7425,self.paraD2, sec = self.FLUT[51])
        h.pt3dadd(-6773.6092,5497.6902,503.9397,self.paraD2, sec = self.FLUT[51])

        h.pt3dadd(-6773.6092,5497.6902,503.9397,self.paraD1, sec = self.MYSA[51])
        h.pt3dadd(-6771.6429,5499.4748,505.3356,self.paraD1, sec = self.MYSA[51])

        h.pt3dadd(-6771.6429,5499.4748,505.3356,self.nodeD, sec = self.node[26])
        h.pt3dadd(-6770.9874,5500.0697,505.8009,self.nodeD, sec = self.node[26])








