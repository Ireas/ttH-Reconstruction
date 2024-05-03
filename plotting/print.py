import sys
import numpy as np
import uproot
from module_table import *


def main():
    # access given root file
    assert len(sys.argv)==2, "Error: root file must be given as only argument"
    root_file = uproot.open(sys.argv[1])
    root_tree = root_file["matched"]

    r_n_jets = root_tree["number_of_jets"].array()
    r_m_t = root_tree["reconstructed_t_m"].array()*1e-3
    r_m_tbar = root_tree["reconstructed_tbar_m"].array()*1e-3
    r_m_w = root_tree["reconstructed_W_from_t_m"].array()*1e-3
    r_m_wbar = root_tree["reconstructed_W_from_tbar_m"].array()*1e-3


    r_en_reco = root_tree["reco_event_number"].array()
    r_en_truth = root_tree["truth_event_number"].array()
    r_cn_reco = root_tree["reco_mc_channel_number"].array()
    r_cn_truth = root_tree["truth_mc_channel_number"].array()

    for (en_reco, en_truth, cn_reco, cn_truth) in zip(r_en_reco, r_en_truth, r_cn_reco, r_cn_truth):
        if(en_reco!=en_truth):
            print("Not same event number in: ", en_reco, " ", en_truth)
            exit()
        if(cn_reco!=cn_truth):
            print("Not same channel number in: ", cn_reco, " ", cn_truth)
            exit()


    total = 0
    succ_t = 0
    succ_tbar = 0
    succ_w = 0
    succ_wbar = 0
    succ_event = 0

    for (n_jets, m_t, m_tbar, m_w, m_wbar) in zip(r_n_jets, r_m_t, r_m_tbar, r_m_w, r_m_wbar):
        total+= 1
        if m_t>1:
            succ_t+= 1
        if m_tbar>1:
            succ_tbar+= 1
        if m_w>1:
            succ_w+= 1
        if m_wbar>1:
            succ_wbar+= 1
        if m_t>1 and m_tbar>1:
            succ_event+= 1
        
    print(min(r_n_jets))
    print(" Success Rates per Object")
    print(" >> Yield: ", total)
    print(" >> event: ", np.round(succ_event/total, 3) )
    print(" >> t:     ", np.round(succ_t/total, 3) )
    print(" >> tbar:  ", np.round(succ_tbar/total, 3) )
    print(" >> w:     ", np.round(succ_w/total, 3) )
    print(" >> wbar:  ", np.round(succ_wbar/total, 3) )



    r_h_dm = root_tree["higgs_decay_mode_custom"].array()
    r_h_d1_d1_pdgids = root_tree["higgs_decay1_decay1_filtered_pdgid"].array()
    r_h_d1_d2_pdgids = root_tree["higgs_decay1_decay2_filtered_pdgid"].array()
    r_h_d2_d1_pdgids = root_tree["higgs_decay2_decay1_filtered_pdgid"].array()
    r_h_d2_d2_pdgids = root_tree["higgs_decay2_decay2_filtered_pdgid"].array()

    total = 0
    succ_t = 0
    succ_tbar = 0
    succ_w = 0
    succ_wbar = 0
    succ_event = 0


	### everything with 8+ jets
    for (n_jets, m_t, m_tbar, m_w, m_wbar, h_dm, pdgid_11, pdgid_12, pdgid_21, pdgid_22) in zip(r_n_jets, r_m_t, r_m_tbar, r_m_w, r_m_wbar, r_h_dm, r_h_d1_d1_pdgids, r_h_d1_d2_pdgids, r_h_d2_d1_pdgids, r_h_d2_d2_pdgids):
        
        if n_jets<8:
            continue

        # full invalid
        if (pdgid_11==-1 or pdgid_12==-1 or pdgid_21==-1 or pdgid_22==-1):
            continue
        
        total+= 1
        
        if m_t>1:
            succ_t+= 1
        if m_tbar>1:
            succ_tbar+= 1
        if m_w>1:
            succ_w+= 1
        if m_wbar>1:
            succ_wbar+= 1
        if m_t>1 and m_tbar>1:
            succ_event+= 1

    print()
    print(" Success Rates per Object for >=8 jets")
    print(" >> Yield: ", total)
    if (total!=0):
        print(" >> event: ", np.round(succ_event/total, 3) )
        print(" >> t:     ", np.round(succ_t/total, 3) )
        print(" >> tbar:  ", np.round(succ_tbar/total, 3) )
        print(" >> w:     ", np.round(succ_w/total, 3) )
        print(" >> wbar:  ", np.round(succ_wbar/total, 3) )

    total = 0
    succ_t = 0
    succ_tbar = 0
    succ_w = 0
    succ_wbar = 0
    succ_event = 0
	
    print("EVERYTHING IS 8+ JETS NOW!")

    for (n_jets, m_t, m_tbar, m_w, m_wbar, h_dm, pdgid_11, pdgid_12, pdgid_21, pdgid_22) in zip(r_n_jets, r_m_t, r_m_tbar, r_m_w, r_m_wbar, r_h_dm, r_h_d1_d1_pdgids, r_h_d1_d2_pdgids, r_h_d2_d1_pdgids, r_h_d2_d2_pdgids):
        if h_dm!=5:
            continue
        if n_jets<8:
            continue
        

        # full invalid
        if (pdgid_11==-1 or pdgid_12==-1 or pdgid_21==-1 or pdgid_22==-1):
            continue
        
        total+= 1
        
        if m_t>1:
            succ_t+= 1
        if m_tbar>1:
            succ_tbar+= 1
        if m_w>1:
            succ_w+= 1
        if m_wbar>1:
            succ_wbar+= 1
        if m_t>1 and m_tbar>1:
            succ_event+= 1

    print()
    print(" Success Rates per Object, H->WW only")
    print(" >> Yield: ", total)
    if (total!=0):
        print(" >> event: ", np.round(succ_event/total, 3) )
        print(" >> t:     ", np.round(succ_t/total, 3) )
        print(" >> tbar:  ", np.round(succ_tbar/total, 3) )
        print(" >> w:     ", np.round(succ_w/total, 3) )
        print(" >> wbar:  ", np.round(succ_wbar/total, 3) )


    total = 0
    succ_t = 0
    succ_tbar = 0
    succ_w = 0
    succ_wbar = 0
    succ_event = 0


    for (n_jets, m_t, m_tbar, m_w, m_wbar, h_dm, pdgid_11, pdgid_12, pdgid_21, pdgid_22) in zip(r_n_jets, r_m_t, r_m_tbar, r_m_w, r_m_wbar, r_h_dm, r_h_d1_d1_pdgids, r_h_d1_d2_pdgids, r_h_d2_d1_pdgids, r_h_d2_d2_pdgids):
        if h_dm!=0:
            continue
        if n_jets<8:
            continue
        
        total+= 1
        
        if m_t>1:
            succ_t+= 1
        if m_tbar>1:
            succ_tbar+= 1
        if m_w>1:
            succ_w+= 1
        if m_wbar>1:
            succ_wbar+= 1
        if m_t>1 and m_tbar>1:
            succ_event+= 1

    print()
    print(" Success Rates per Object, H-bb only")
    print(" >> Yield: ", total)
    if (total!=0):
        print(" >> event: ", np.round(succ_event/total, 3) )
        print(" >> t:     ", np.round(succ_t/total, 3) )
        print(" >> tbar:  ", np.round(succ_tbar/total, 3) )
        print(" >> w:     ", np.round(succ_w/total, 3) )
        print(" >> wbar:  ", np.round(succ_wbar/total, 3) )

    total = 0
    succ_t = 0
    succ_tbar = 0
    succ_w = 0
    succ_wbar = 0
    succ_event = 0

    for (n_jets, m_t, m_tbar, m_w, m_wbar, h_dm, pdgid_11, pdgid_12, pdgid_21, pdgid_22) in zip(r_n_jets, r_m_t, r_m_tbar, r_m_w, r_m_wbar, r_h_dm, r_h_d1_d1_pdgids, r_h_d1_d2_pdgids, r_h_d2_d1_pdgids, r_h_d2_d2_pdgids):
        # ww-decay
        if h_dm!=5:
            continue
        if n_jets<8:
            continue

        # full invalid
        if (pdgid_11==-1 or pdgid_12==-1 or pdgid_21==-1 or pdgid_22==-1):
            continue
        
        # full hadronic
        if not ( (pdgid_11>0 and pdgid_11<9) and (pdgid_12>0 and pdgid_12<9) and (pdgid_21>0 and pdgid_21<9) and (pdgid_22>0 and pdgid_22<9) ):
            continue


        total+= 1
        if m_t>1:
            succ_t+= 1
        if m_tbar>1:
            succ_tbar+= 1
        if m_w>1:
            succ_w+= 1
        if m_wbar>1:
            succ_wbar+= 1
        if m_t>1 and m_tbar>1:
            succ_event+= 1
        
    print()
    print(" Success Rates per Object, H->WW full hadronic only")
    print(" >> Yield: ", total)
    if (total!=0):
        print(" >> event: ", np.round(succ_event/total, 3) )
        print(" >> t:     ", np.round(succ_t/total, 3) )
        print(" >> tbar:  ", np.round(succ_tbar/total, 3) )
        print(" >> w:     ", np.round(succ_w/total, 3) )
        print(" >> wbar:  ", np.round(succ_wbar/total, 3) )


    total = 0
    succ_t = 0
    succ_tbar = 0
    succ_w = 0
    succ_wbar = 0
    succ_event = 0

    for (n_jets, m_t, m_tbar, m_w, m_wbar, h_dm, pdgid_11, pdgid_12, pdgid_21, pdgid_22) in zip(r_n_jets, r_m_t, r_m_tbar, r_m_w, r_m_wbar, r_h_dm, r_h_d1_d1_pdgids, r_h_d1_d2_pdgids, r_h_d2_d1_pdgids, r_h_d2_d2_pdgids):
        # ww-decay
        if h_dm!=5:
            continue
        if n_jets<8:
            continue

        # full invalid
        if (pdgid_11==-1 or pdgid_12==-1 or pdgid_21==-1 or pdgid_22==-1):
            continue
        
        # full hadronic
        if ( (pdgid_11>0 and pdgid_11<9) and (pdgid_12>0 and pdgid_12<9) and (pdgid_21>0 and pdgid_21<9) and (pdgid_22>0 and pdgid_22<9) ):
            continue

        # full leptonic
        if ( (pdgid_11>10 and pdgid_11<19) and (pdgid_12>10 and pdgid_12<19) and (pdgid_21>10 and pdgid_21<19) and (pdgid_22>10 and pdgid_22<19) ):
            continue


        # implement check if qq / ev come from same w-boson 

        total+= 1
        if m_t>1:
            succ_t+= 1
        if m_tbar>1:
            succ_tbar+= 1
        if m_w>1:
            succ_w+= 1
        if m_wbar>1:
            succ_wbar+= 1
        if m_t>1 and m_tbar>1:
            succ_event+= 1
        
    print()
    print(" Success Rates per Object, H->WW semi-leptonic only")
    print(" >> Yield: ", total)
    if (total!=0):
        print(" >> event: ", np.round(succ_event/total, 3) )
        print(" >> t:     ", np.round(succ_t/total, 3) )
        print(" >> tbar:  ", np.round(succ_tbar/total, 3) )
        print(" >> w:     ", np.round(succ_w/total, 3) )
        print(" >> wbar:  ", np.round(succ_wbar/total, 3) )

    total = 0
    succ_t = 0
    succ_tbar = 0
    succ_w = 0
    succ_wbar = 0
    succ_event = 0

    for (n_jets, m_t, m_tbar, m_w, m_wbar, h_dm, pdgid_11, pdgid_12, pdgid_21, pdgid_22) in zip(r_n_jets, r_m_t, r_m_tbar, r_m_w, r_m_wbar, r_h_dm, r_h_d1_d1_pdgids, r_h_d1_d2_pdgids, r_h_d2_d1_pdgids, r_h_d2_d2_pdgids):
        # ww-decay
        if h_dm!=5:
            continue
        if n_jets<8:
            continue

        # full invalid
        if (pdgid_11==-1 or pdgid_12==-1 or pdgid_21==-1 or pdgid_22==-1):
            continue

        # full leptonic
        if not ( (pdgid_11>10 and pdgid_11<19) and (pdgid_12>10 and pdgid_12<19) and (pdgid_21>10 and pdgid_21<19) and (pdgid_22>10 and pdgid_22<19) ):
            continue

        total+= 1
        if m_t>1:
            succ_t+= 1
        if m_tbar>1:
            succ_tbar+= 1
        if m_w>1:
            succ_w+= 1
        if m_wbar>1:
            succ_wbar+= 1
        if m_t>1 and m_tbar>1:
            succ_event+= 1

        
    print()
    print(" Success Rates per Object, H->WW full-leptonic only")
    print(" >> Yield: ", total)
    if (total!=0):
        print(" >> event: ", np.round(succ_event/total, 3) )
        print(" >> t:     ", np.round(succ_t/total, 3) )
        print(" >> tbar:  ", np.round(succ_tbar/total, 3) )
        print(" >> w:     ", np.round(succ_w/total, 3) )
        print(" >> wbar:  ", np.round(succ_wbar/total, 3) )


if __name__=='__main__':
	main()	
