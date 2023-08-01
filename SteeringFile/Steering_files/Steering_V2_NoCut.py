import basf2 as b2
import modularAnalysis as ma
from variables import variables as vm
import variables.collections as vc
import variables.utils as vu
import kinfit as vx
import vertex as verf
from stdPhotons import stdPhotons
import os, sys
import pdg

outfile=sys.argv[1]
runtype=sys.argv[2]

mypath = b2.create_path()
pdg.add_particle('ap',53,0,0,0,0)

ma.inputMdstList(environmentType='default', filelist=runtype, path=mypath)

ma.fillParticleList('mu+:select', '[abs(dz) < 5.0] and [abs(dr) < 2.0] and pt > 0.2', path=mypath)
ma.fillParticleList('gamma:select', '[useCMSFrame(E) > 0.2 and thetaInCDCAcceptance]', path=mypath)
apcut = ''

ma.reconstructDecay('ap:mumu -> mu-:select mu+:select',apcut,path =mypath)

#ma.applyCuts('ap:mumu','[daughter(0,muonID) > 0.6] and [daughter(1,muonID) > 0.6]',path=mypath)
#ma.applyCuts('ap:mumu','[daughter(0,nPXDHits) > 0] and [daughter(1,nPXDHits) > 0]',path=mypath)
#ma.applyCuts('ap:mumu','theta > 0.175',path=mypath)
#ma.applyCuts('mu-:select','cosTheta > -0.5',path=mypath)

vm.addAlias('cms_cosTheta','useCMSFrame(cosTheta)')

prefitvars = ['px', 'py', 'pz', 'pt', 'p', 'E','cosTheta','cms_cosTheta','dz','dr', 'dphi','M']
prefitlist = [var+'_pf' for var in prefitvars]
prefitdict = dict(zip(prefitvars,prefitlist))
ma.variablesToExtraInfo('ap:mumu', variables = prefitdict, path=mypath)
for item in prefitlist:
    vm.addAlias(item, 'extraInfo('+item+')')
    print(item)

verf.kFit('ap:mumu', 0.0 , path=mypath)

ma.reconstructDecay('vpho:Total -> ap:mumu gamma:select','',path =mypath)
vx.fitKinematic4C('vpho:Total',daughtersUpdate=False ,path=mypath)
ma.matchMCTruth('vpho:Total', mypath)

#########Variables##########

vm.addAlias('cms_p','useCMSFrame(p)')
vm.addAlias('cms_theta','useCMSFrame(theta)')
vm.addAlias('cms_E','useCMSFrame(E)')
vm.addAlias('deltaM','massDifference(0)')
vm.addAlias('deltaMErr','massDifferenceError(0)')
vm.addAlias('dau_theta','daughterAngle(0, 1)')
vm.addAlias('dau_phi','daughterDiffOf(0, 1, phi)')
vm.addAlias('dau_hel1','cosHelicityAngle(0,0)')
vm.addAlias('dau_hel2','cosHelicityAngle(0,1)')

cmsvars = ['cms_p','cms_E','cms_cosTheta','cms_theta']
basevars = ['theta','cosTheta','phi']+vc.inv_mass + vc.kinematics + vc.vertex + vc.track + vc.mc_truth + cmsvars + vc.trackfit_parameters
  
vphovars = basevars + ['mcPDG','charge','dau_hel1','dau_hel2']
apvars = basevars +['dau_theta','dau_phi','pValue','d0','phi0','z0','vertexDistance','flightDistance']
phovars = basevars+['clusterErrorTiming','clusterE1E9','inCDCAcceptance']
dauvars = basevars + ['electronID','muonID','pionID','nPXDHits','d0', 'z0', 'pValue', 'charge']


var_lst  = vu.create_aliases(vphovars,'{variable}','vpho')
var_lst += vu.create_aliases(apvars+prefitlist,  'daughter(0,{variable})','ap')
var_lst += vu.create_aliases(phovars, 'daughter(1,{variable})','pho')
var_lst += vu.create_aliases(dauvars, 'daughter(0,daughter(0,{variable}))','d1')
var_lst += vu.create_aliases(dauvars, 'daughter(0,daughter(1,{variable}))','d2')

#########NTUPLE##########

ma.variablesToNtuple('vpho:Total', var_lst , filename=outfile, treename='mumu' ,path=mypath)

b2.process(path = mypath)



