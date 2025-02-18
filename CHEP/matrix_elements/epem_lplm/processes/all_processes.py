from __future__ import division
from model.aloha_methods import *
from model.wavefunctions import *
class Matrix_2_epem_epem_no_z(object):

    def __init__(self):
        """define the object"""
        self.clean()

    def clean(self):
        self.jamp = []

    def get_external_masses(self, model):

        return ( (model.ZERO, model.ZERO), (model.ZERO, model.ZERO) )

    def smatrix(self,p, model):
        #  
        #  MadGraph5_aMC@NLO v. 3.5.7, 2024-11-29
        #  By the MadGraph5_aMC@NLO Development Team
        #  Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
        # 
        # MadGraph5_aMC@NLO StandAlone Version
        # 
        # Returns amplitude squared summed/avg over colors
        # and helicities
        # for the point in phase space P(0:3,NEXTERNAL)
        #  
        # Process: e+ e- > e+ e- WEIGHTED<=4 / z @2
        #  
        # Clean additional output
        #
        self.clean()
        #  
        # CONSTANTS
        #  
        nexternal = 4
        ndiags = 2
        ncomb = 16
        #  
        # LOCAL VARIABLES 
        #  
        helicities = [ \
        [-1,1,1,-1],
        [-1,1,1,1],
        [-1,1,-1,-1],
        [-1,1,-1,1],
        [-1,-1,1,-1],
        [-1,-1,1,1],
        [-1,-1,-1,-1],
        [-1,-1,-1,1],
        [1,1,1,-1],
        [1,1,1,1],
        [1,1,-1,-1],
        [1,1,-1,1],
        [1,-1,1,-1],
        [1,-1,1,1],
        [1,-1,-1,-1],
        [1,-1,-1,1]]
        denominator = 4
        # ----------
        # BEGIN CODE
        # ----------
        self.amp2 = [0.] * ndiags
        self.helEvals = []
        ans = 0.
        for hel in helicities:
            t = self.matrix(p, hel, model)
            ans = ans + t
            self.helEvals.append([hel, t.real / denominator ])
        ans = ans / denominator
        return ans.real

    def matrix(self, p, hel, model):
        #  
        #  MadGraph5_aMC@NLO v. 3.5.7, 2024-11-29
        #  By the MadGraph5_aMC@NLO Development Team
        #  Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
        #
        # Returns amplitude squared summed/avg over colors
        # for the point with external lines W(0:6,NEXTERNAL)
        #
        # Process: e+ e- > e+ e- WEIGHTED<=4 / z @2
        #  
        #  
        # Process parameters
        #  
        ngraphs = 2
        nexternal = 4
        nwavefuncs = 5
        ncolor = 1
        ZERO = 0.
        #  
        # Color matrix
        #  
        denom = [1];
        cf = [[1]];
        #
        # Model parameters
        #
        
        GC_3 = model.GC_3
        # ----------
        # Begin code
        # ----------
        amp = [None] * ngraphs
        w = [None] * nwavefuncs
        w[0] = oxxxxx(p[0],ZERO,hel[0],-1)
        w[1] = ixxxxx(p[1],ZERO,hel[1],+1)
        w[2] = ixxxxx(p[2],ZERO,hel[2],-1)
        w[3] = oxxxxx(p[3],ZERO,hel[3],+1)
        w[4]= FFV1P0_3(w[1],w[0],GC_3,ZERO,ZERO)
        # Amplitude(s) for diagram number 1
        amp[0]= FFV1_0(w[2],w[3],w[4],GC_3)
        w[4]= FFV1P0_3(w[2],w[0],GC_3,ZERO,ZERO)
        # Amplitude(s) for diagram number 2
        amp[1]= FFV1_0(w[1],w[3],w[4],GC_3)

        jamp = [None] * ncolor

        jamp[0] = -amp[0]+amp[1]

        self.amp2[0]+=abs(amp[0]*amp[0].conjugate())
        self.amp2[1]+=abs(amp[1]*amp[1].conjugate())
        matrix = 0.
        for i in range(ncolor):
            ztemp = 0
            for j in range(ncolor):
                ztemp = ztemp + cf[i][j]*jamp[j]
            matrix = matrix + ztemp * jamp[i].conjugate()/denom[i]   
        self.jamp.append(jamp)

        return matrix

class Matrix_1_epem_mupmum_no_z(object):

    def __init__(self):
        """define the object"""
        self.clean()

    def clean(self):
        self.jamp = []

    def get_external_masses(self, model):

        return ( (model.ZERO, model.ZERO), (model.ZERO, model.ZERO) )

    def smatrix(self,p, model):
        #  
        #  MadGraph5_aMC@NLO v. 3.5.7, 2024-11-29
        #  By the MadGraph5_aMC@NLO Development Team
        #  Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
        # 
        # MadGraph5_aMC@NLO StandAlone Version
        # 
        # Returns amplitude squared summed/avg over colors
        # and helicities
        # for the point in phase space P(0:3,NEXTERNAL)
        #  
        # Process: e+ e- > mu+ mu- WEIGHTED<=4 / z @1
        #  
        # Clean additional output
        #
        self.clean()
        #  
        # CONSTANTS
        #  
        nexternal = 4
        ndiags = 1
        ncomb = 16
        #  
        # LOCAL VARIABLES 
        #  
        helicities = [ \
        [-1,1,1,-1],
        [-1,1,1,1],
        [-1,1,-1,-1],
        [-1,1,-1,1],
        [-1,-1,1,-1],
        [-1,-1,1,1],
        [-1,-1,-1,-1],
        [-1,-1,-1,1],
        [1,1,1,-1],
        [1,1,1,1],
        [1,1,-1,-1],
        [1,1,-1,1],
        [1,-1,1,-1],
        [1,-1,1,1],
        [1,-1,-1,-1],
        [1,-1,-1,1]]
        denominator = 4
        # ----------
        # BEGIN CODE
        # ----------
        self.amp2 = [0.] * ndiags
        self.helEvals = []
        ans = 0.
        for hel in helicities:
            t = self.matrix(p, hel, model)
            ans = ans + t
            self.helEvals.append([hel, t.real / denominator ])
        ans = ans / denominator
        return ans.real

    def matrix(self, p, hel, model):
        #  
        #  MadGraph5_aMC@NLO v. 3.5.7, 2024-11-29
        #  By the MadGraph5_aMC@NLO Development Team
        #  Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
        #
        # Returns amplitude squared summed/avg over colors
        # for the point with external lines W(0:6,NEXTERNAL)
        #
        # Process: e+ e- > mu+ mu- WEIGHTED<=4 / z @1
        #  
        #  
        # Process parameters
        #  
        ngraphs = 1
        nexternal = 4
        nwavefuncs = 5
        ncolor = 1
        ZERO = 0.
        #  
        # Color matrix
        #  
        denom = [1];
        cf = [[1]];
        #
        # Model parameters
        #
        
        GC_3 = model.GC_3
        # ----------
        # Begin code
        # ----------
        amp = [None] * ngraphs
        w = [None] * nwavefuncs
        w[0] = oxxxxx(p[0],ZERO,hel[0],-1)
        w[1] = ixxxxx(p[1],ZERO,hel[1],+1)
        w[2] = ixxxxx(p[2],ZERO,hel[2],-1)
        w[3] = oxxxxx(p[3],ZERO,hel[3],+1)
        w[4]= FFV1P0_3(w[1],w[0],GC_3,ZERO,ZERO)
        # Amplitude(s) for diagram number 1
        amp[0]= FFV1_0(w[2],w[3],w[4],GC_3)

        jamp = [None] * ncolor

        jamp[0] = -amp[0]

        self.amp2[0]+=abs(amp[0]*amp[0].conjugate())
        matrix = 0.
        for i in range(ncolor):
            ztemp = 0
            for j in range(ncolor):
                ztemp = ztemp + cf[i][j]*jamp[j]
            matrix = matrix + ztemp * jamp[i].conjugate()/denom[i]   
        self.jamp.append(jamp)

        return matrix

