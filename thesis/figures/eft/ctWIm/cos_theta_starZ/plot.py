####  Header  #####
####################
import ROOT as root
from ROOT import TFile
import atlasplots as aplt



#####  Setup Style  ######
##########################
aplt.set_atlas_style()

### Style
width = 800
height = 800

offsetA = 5.0
offsetB = 3.0
offsetC = 2.0
ratiolimA = 0.30
ratiolimB = 0.15
ratiolimC = 0.10

string_xlabel = r"cos(\theta*) for the \it{Z}-boson"
string_file = "cos_theta_starZ"

### Create Figures
figA, (ax1A, ax2A) = aplt.ratio_plot(name="figA_ttZ", figsize=(width, height), hspace=0.05)
figB, (ax1B, ax2B) = aplt.ratio_plot(name="figB_ttW", figsize=(width, height), hspace=0.05)
figC, (ax1C, ax2C) = aplt.ratio_plot(name="figC_Fakes", figsize=(width, height), hspace=0.05)



#####  Load Histograms  #####
#############################
smeft_file = TFile("output.root")

### Histograms
ttZ_default = smeft_file.Get("ttZ_default")
ttZ_c1 = smeft_file.Get("ttZ_c1")
ttZ_c2 = smeft_file.Get("ttZ_c2")
ttZ_c3 = smeft_file.Get("ttZ_c3")
ttZ_c4 = smeft_file.Get("ttZ_c4")

ttW_default = smeft_file.Get("ttW_default")
ttW_c1 = smeft_file.Get("ttW_c1")
ttW_c2 = smeft_file.Get("ttW_c2")
ttW_c3 = smeft_file.Get("ttW_c3")
ttW_c4 = smeft_file.Get("ttW_c4")

Fakes_default = smeft_file.Get("Fakes_default")
Fakes_c1 = smeft_file.Get("Fakes_c1")
Fakes_c2 = smeft_file.Get("Fakes_c2")
Fakes_c3 = smeft_file.Get("Fakes_c3")
Fakes_c4 = smeft_file.Get("Fakes_c4")



#####  Define Objects Main  #####
#################################

### Main Cuts
ylim_A = max( aplt.root_helpers.hist_max(ttZ_default), aplt.root_helpers.hist_max(ttZ_c1), aplt.root_helpers.hist_max(ttZ_c2), aplt.root_helpers.hist_max(ttZ_c3), aplt.root_helpers.hist_max(ttZ_c4))
ylim_B = max( aplt.root_helpers.hist_max(ttW_default), aplt.root_helpers.hist_max(ttW_c1), aplt.root_helpers.hist_max(ttW_c2), aplt.root_helpers.hist_max(ttW_c3), aplt.root_helpers.hist_max(ttW_c4) )
ylim_C = max( aplt.root_helpers.hist_max(Fakes_default), aplt.root_helpers.hist_max(Fakes_c1), aplt.root_helpers.hist_max(Fakes_c2), aplt.root_helpers.hist_max(Fakes_c3), aplt.root_helpers.hist_max(Fakes_c4) )


#####  Plot ax1 (Main)  #####
#############################
ax1A.plot(ttZ_default, "E0", markerstyle=1, linewidth=2, linecolor=root.kBlack)
ax1A.plot(ttZ_c1, "E0", markerstyle=1, linewidth=2, linecolor=root.kRed-4)
ax1A.plot(ttZ_c2, "E0", markerstyle=1, linewidth=2, linecolor=root.kBlue-4)
ax1A.plot(ttZ_c3, "E0", markerstyle=1, linewidth=2, linecolor=root.kGreen)
ax1A.plot(ttZ_c4, "E0", markerstyle=1, linewidth=2, linecolor=root.kMagenta)
err_ttZ_c1 = aplt.root_helpers.hist_to_graph(ttZ_c1, show_bin_width=True)
err_ttZ_c2 = aplt.root_helpers.hist_to_graph(ttZ_c2, show_bin_width=True)
err_ttZ_c3 = aplt.root_helpers.hist_to_graph(ttZ_c3, show_bin_width=True)
err_ttZ_c4 = aplt.root_helpers.hist_to_graph(ttZ_c4, show_bin_width=True)
ax1A.plot(err_ttZ_c1, "2", fillcolor=root.kRed-4, fillstyle=3254, linewidth=0)
ax1A.plot(err_ttZ_c2, "2", fillcolor=root.kBlue-4, fillstyle=3245, linewidth=0)
ax1A.plot(err_ttZ_c3, "2", fillcolor=root.kGreen, fillstyle=3544, linewidth=0)
ax1A.plot(err_ttZ_c4, "2", fillcolor=root.kMagenta, fillstyle=3359, linewidth=0)


ax1B.plot(ttW_default, "E0", markerstyle=1, linewidth=2, linecolor=root.kBlack )
ax1B.plot(ttW_c1, "E0", markerstyle=1, linewidth=2, linecolor=root.kRed-4)
ax1B.plot(ttW_c2, "E0", markerstyle=1, linewidth=2, linecolor=root.kBlue-4)
ax1B.plot(ttW_c3, "E0", markerstyle=1, linewidth=2, linecolor=root.kGreen)
ax1B.plot(ttW_c4, "E0", markerstyle=1, linewidth=2, linecolor=root.kMagenta)
err_ttW_c1 = aplt.root_helpers.hist_to_graph(ttW_c1, show_bin_width=True)
err_ttW_c2 = aplt.root_helpers.hist_to_graph(ttW_c2, show_bin_width=True)
err_ttW_c3 = aplt.root_helpers.hist_to_graph(ttW_c3, show_bin_width=True)
err_ttW_c4 = aplt.root_helpers.hist_to_graph(ttW_c4, show_bin_width=True)
ax1B.plot(err_ttW_c1, "2", fillcolor=root.kRed-4, fillstyle=3254, linewidth=0)
ax1B.plot(err_ttW_c2, "2", fillcolor=root.kBlue-4, fillstyle=3245, linewidth=0)
ax1B.plot(err_ttW_c3, "2", fillcolor=root.kGreen, fillstyle=3544, linewidth=0)
ax1B.plot(err_ttW_c4, "2", fillcolor=root.kMagenta, fillstyle=3359, linewidth=0)


ax1C.plot(Fakes_default, "E0", markerstyle=1, linewidth=2, linecolor=root.kBlack )
ax1C.plot(Fakes_c1, "E0", markerstyle=1, linewidth=2, linecolor=root.kRed-4)
ax1C.plot(Fakes_c2, "E0", markerstyle=1, linewidth=2, linecolor=root.kBlue-4)
ax1C.plot(Fakes_c3, "E0", markerstyle=1, linewidth=2, linecolor=root.kGreen)
ax1C.plot(Fakes_c4, "E0", markerstyle=1, linewidth=2, linecolor=root.kMagenta)
err_Fakes_c1 = aplt.root_helpers.hist_to_graph(Fakes_c1, show_bin_width=True)
err_Fakes_c2 = aplt.root_helpers.hist_to_graph(Fakes_c2, show_bin_width=True)
err_Fakes_c3 = aplt.root_helpers.hist_to_graph(Fakes_c3, show_bin_width=True)
err_Fakes_c4 = aplt.root_helpers.hist_to_graph(Fakes_c4, show_bin_width=True)
ax1C.plot(err_Fakes_c1, "2", fillcolor=root.kRed-4, fillstyle=3254, linewidth=0)
ax1C.plot(err_Fakes_c2, "2", fillcolor=root.kBlue-4, fillstyle=3245, linewidth=0)
ax1C.plot(err_Fakes_c3, "2", fillcolor=root.kGreen, fillstyle=3544, linewidth=0)
ax1C.plot(err_Fakes_c4, "2", fillcolor=root.kMagenta, fillstyle=3359, linewidth=0)


ax1A.set_ylim(0, ylim_A+offsetA)
ax1B.set_ylim(0, ylim_B+offsetB)
ax1C.set_ylim(0, ylim_C+offsetC)
#####  Define Objects Ratio  #####
##################################

### Horizontal Ratio Lines
ax2A.set_xlim(ax1A.get_xlim())
ax2B.set_xlim(ax1B.get_xlim())
ax2C.set_xlim(ax1C.get_xlim())

lineA = root.TLine(ax1A.get_xlim()[0], 1, ax1A.get_xlim()[1], 1)
lineB = root.TLine(ax1B.get_xlim()[0], 1, ax1B.get_xlim()[1], 1)
lineC = root.TLine(ax1C.get_xlim()[0], 1, ax1C.get_xlim()[1], 1)

ttZ_default_clone = ttZ_default.Clone("ttZ_default_clone")
ratioAc1 = ttZ_c1.Clone("ratioAc1")
ratioAc2 = ttZ_c2.Clone("ratioAc2")
ratioAc3 = ttZ_c3.Clone("ratioAc3")
ratioAc4 = ttZ_c4.Clone("ratioAc4")
ratioAc1.Divide(ttZ_default_clone)
ratioAc2.Divide(ttZ_default_clone)
ratioAc3.Divide(ttZ_default_clone)
ratioAc4.Divide(ttZ_default_clone)
ratioAc1_graph = aplt.root_helpers.hist_to_graph(ratioAc1, show_bin_width=True)
ratioAc2_graph = aplt.root_helpers.hist_to_graph(ratioAc2, show_bin_width=True)
ratioAc3_graph = aplt.root_helpers.hist_to_graph(ratioAc3, show_bin_width=True)
ratioAc4_graph = aplt.root_helpers.hist_to_graph(ratioAc4, show_bin_width=True)

ttW_default_clone = ttW_default.Clone("ttW_default_clone")
ratioBc1 = ttW_c1.Clone("ratioBc1")
ratioBc2 = ttW_c2.Clone("ratioBc2")
ratioBc3 = ttW_c3.Clone("ratioBc3")
ratioBc4 = ttW_c4.Clone("ratioBc4")
ratioBc1.Divide(ttW_default_clone)
ratioBc2.Divide(ttW_default_clone)
ratioBc3.Divide(ttW_default_clone)
ratioBc4.Divide(ttW_default_clone)
ratioBc1_graph = aplt.root_helpers.hist_to_graph(ratioBc1, show_bin_width=True)
ratioBc2_graph = aplt.root_helpers.hist_to_graph(ratioBc2, show_bin_width=True)
ratioBc3_graph = aplt.root_helpers.hist_to_graph(ratioBc3, show_bin_width=True)
ratioBc4_graph = aplt.root_helpers.hist_to_graph(ratioBc4, show_bin_width=True)

Fakes_default_clone = Fakes_default.Clone("Fakes_default_clone")
ratioCc1 = Fakes_c1.Clone("ratioCc1")
ratioCc2 = Fakes_c2.Clone("ratioCc2")
ratioCc3 = Fakes_c3.Clone("ratioCc3")
ratioCc4 = Fakes_c4.Clone("ratioCc4")
ratioCc1.Divide(Fakes_default_clone)
ratioCc2.Divide(Fakes_default_clone)
ratioCc3.Divide(Fakes_default_clone)
ratioCc4.Divide(Fakes_default_clone)
ratioCc1_graph = aplt.root_helpers.hist_to_graph(ratioCc1, show_bin_width=True)
ratioCc2_graph = aplt.root_helpers.hist_to_graph(ratioCc2, show_bin_width=True)
ratioCc3_graph = aplt.root_helpers.hist_to_graph(ratioCc3, show_bin_width=True)
ratioCc4_graph = aplt.root_helpers.hist_to_graph(ratioCc4, show_bin_width=True)


#####  Plot ax2 (Ratio)  #####
##############################
ax2A.plot(lineA)
ax2B.plot(lineB)
ax2C.plot(lineC)

ax2A.plot(ratioAc1, "0", linewidth=2, linecolor=root.kRed-4)
ax2A.plot(ratioAc2, "0", linewidth=2, linecolor=root.kBlue-4)
ax2A.plot(ratioAc3, "0", linewidth=2, linecolor=root.kGreen)
ax2A.plot(ratioAc4, "0", linewidth=2, linecolor=root.kMagenta)
ax2A.plot(ratioAc1_graph, "2", fillcolor=root.kRed-4, fillstyle=3254, linewidth=0)
ax2A.plot(ratioAc2_graph, "2", fillcolor=root.kBlue-4, fillstyle=3245, linewidth=0)
ax2A.plot(ratioAc3_graph, "2", fillcolor=root.kGreen, fillstyle=3544, linewidth=0)
ax2A.plot(ratioAc4_graph, "2", fillcolor=root.kMagenta, fillstyle=3359, linewidth=0)


ax2B.plot(ratioBc1, "0", linewidth=2, linecolor=root.kRed-4)
ax2B.plot(ratioBc2, "0", linewidth=2, linecolor=root.kBlue-4)
ax2B.plot(ratioBc3, "0", linewidth=2, linecolor=root.kGreen)
ax2B.plot(ratioBc4, "0", linewidth=2, linecolor=root.kMagenta)
ax2B.plot(ratioBc1_graph, "2", fillcolor=root.kRed-4, fillstyle=3254, linewidth=0)
ax2B.plot(ratioBc2_graph, "2", fillcolor=root.kBlue-4, fillstyle=3245, linewidth=0)
ax2B.plot(ratioBc3_graph, "2", fillcolor=root.kGreen, fillstyle=3544, linewidth=0)
ax2B.plot(ratioBc4_graph, "2", fillcolor=root.kMagenta, fillstyle=3359, linewidth=0)


ax2C.plot(ratioCc1, "0", linewidth=2, linecolor=root.kRed-4)
ax2C.plot(ratioCc2, "0", linewidth=2, linecolor=root.kBlue-4)
ax2C.plot(ratioCc3, "0", linewidth=2, linecolor=root.kGreen)
ax2C.plot(ratioCc4, "0", linewidth=2, linecolor=root.kMagenta)
ax2C.plot(ratioCc1_graph, "2", fillcolor=root.kRed-4, fillstyle=3254, linewidth=0)
ax2C.plot(ratioCc2_graph, "2", fillcolor=root.kBlue-4, fillstyle=3245, linewidth=0)
ax2C.plot(ratioCc3_graph, "2", fillcolor=root.kGreen, fillstyle=3544, linewidth=0)
ax2C.plot(ratioCc4_graph, "2", fillcolor=root.kMagenta, fillstyle=3359, linewidth=0)



#####  Style Figure  #####
##########################

### Margins
ax1A.add_margins(top=0.25)
ax1B.add_margins(top=0.25)
ax1C.add_margins(top=0.25)

### Axis
ax2A.set_xlabel(string_xlabel)
ax2B.set_xlabel(string_xlabel)
ax2C.set_xlabel(string_xlabel)

ax1A.set_ylabel("Events")
ax1B.set_ylabel("Events")
ax1C.set_ylabel("Events")

ax2A.set_ylabel("Ratio", loc="centre")
ax2B.set_ylabel("Ratio", loc="centre")
ax2B.set_ylabel("Ratio", loc="centre")

ax2A.set_ylim(1.0-ratiolimA, 1+ratiolimA)
ax2B.set_ylim(1.0-ratiolimB, 1+ratiolimB)
ax2C.set_ylim(1.0-ratiolimC, 1+ratiolimC)

#####  Labels  #####
####################

### Labels Plot A
ax1A.cd()
legend = ax1A.legend(loc=(0.7, 0.6, 1-root.gPad.GetRightMargin(), 1-root.gPad.GetTopMargin()-0.05), textsize=22, labelspacing=250 )
legend.AddEntry(ttZ_default,    "Default", "L")
legend.AddEntry(ttZ_c1,         "ctWIm +0.6", "L")
legend.AddEntry(ttZ_c2,         "ctWIm -0.8", "L")
legend.AddEntry(ttZ_c3,         "ctWIm +1.2", "L")
legend.AddEntry(ttZ_c4,         "ctWIm -1.4", "L")
#aplt.atlas_label(text="Internal", loc="upper left")
ax1A.text(0.2, 0.88, "#sqrt{s}=13 TeV, 139 fb^{-1}", size=22, align=13)
ax1A.text(0.2, 0.82, "SR-ttZ-NN", size=22, align=13)
ax1A.text(0.2, 0.76, "Variations on imaginary part of C_{tW}", size=22, align=13)

### Labels Plot B
ax1B.cd()
legend = ax1B.legend(loc=(0.7, 0.6, 1-root.gPad.GetRightMargin(), 1-root.gPad.GetTopMargin()-0.05), textsize=22, labelspacing=250 )
legend.AddEntry(ttW_default,    "Default", "L")
legend.AddEntry(ttW_c1,         "ctWIm +0.6", "L")
legend.AddEntry(ttW_c2,         "ctWIm -0.8", "L")
legend.AddEntry(ttW_c3,         "ctWIm +1.2", "L")
legend.AddEntry(ttW_c4,         "ctWIm -1.4", "L")
#aplt.atlas_label(text="Internal", loc="upper left")
ax1B.text(0.2, 0.88, "#sqrt{s}=13 TeV, 139 fb^{-1}", size=22, align=13)
ax1B.text(0.2, 0.82, "CR-ttW-NN", size=22, align=13)
ax1B.text(0.2, 0.86, "Variations on imaginary part of C_{tW}", size=22, align=13)

### Labels Plot C
ax1C.cd()
legend = ax1C.legend(loc=(0.7, 0.6, 1-root.gPad.GetRightMargin(), 1-root.gPad.GetTopMargin()-0.05), textsize=22, labelspacing=250 )
legend.AddEntry(Fakes_default,    "Default", "L")
legend.AddEntry(Fakes_c1,         "ctWIm +0.6", "L")
legend.AddEntry(Fakes_c2,         "ctWIm -0.8", "L")
legend.AddEntry(Fakes_c3,         "ctWIm +1.2", "L")
legend.AddEntry(Fakes_c4,         "ctWIm -1.4", "L")
#aplt.atlas_label(text="Internal", loc="upper left")
ax1C.text(0.2, 0.88, "#sqrt{s}=13 TeV, 139 fb^{-1}", size=22, align=13)
ax1C.text(0.2, 0.82, "CR-Fakes-NN", size=22, align=13)
ax1C.text(0.2, 0.76, "Variations on imaginary part of C_{tW}", size=22, align=13)



#####  Print Figure  #####
##########################
figA.savefig("_"+string_file+"_region00_ttZ.pdf")
figB.savefig("_"+string_file+"_region01_ttW.pdf")
figC.savefig("_"+string_file+"_region02_Fakes.pdf")
