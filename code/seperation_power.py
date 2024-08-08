import os
import sys
import uproot
import numpy as np
import matplotlib.pyplot as plt


# CONSTANTS
# files
ROOT_INPUT_FILE = "/media/ireas/Data/v4/injected/injection.root"

SHOW_PLOTS = False

# binning
NUMBER_OF_BINS = 20

DISTRIBUTION_VARIABLE_NAMES = [
    "spanet/t1_marginal_probability",
    "spanet/t1_detection_probability",
    "spanet/t1_assignment_probability",
    "spanet/t2_marginal_probability",
    "spanet/t2_detection_probability",
    "spanet/t2_assignment_probability",
    "spanet/HW_marginal_probability",
    "spanet/HW_detection_probability",
    "spanet/HW_assignment_probability",
]

# classifier for signal and background
CLASSIFIER_VARIABLE_NAME = "matched/classifier_event_status"


def main():

    for variable in DISTRIBUTION_VARIABLE_NAMES:
        # get signal and background
        signal, background = get_distributions(variable)

        if (len(signal)==0 or len(background)==0 ):
            print(">> signal", len(signal), signal)
            print(">> background", len(background), background)
            print(">> empty entries, skipping")
            continue
        
        # print separation power value
        print( calculate_separation_power(signal, background) )

        # show separation plot
        plot_separation(signal, background, variable)
    



def get_distributions(variable):
    signal = np.array([])
    background = np.array([])

    # access root file
    print(f"> {variable}")
    with uproot.open(ROOT_INPUT_FILE) as root_file:			
        classifier = np.array( root_file[CLASSIFIER_VARIABLE_NAME].array() )
        distribution = np.array( root_file[variable].array() )

        # filter by event type
        for (classification, value) in zip(classifier, distribution):
            if classification==-1:
                background = np.append(background, [value])
            else:
                signal = np.append(signal, [value])

    return signal, background


def transform_to_histogram(signal, background, normalize=False):
    # get binning parameters
    min_value = np.min([np.min(signal), np.min(background)])
    max_value = np.max([np.max(signal), np.max(background)])
    bins = np.linspace(min_value, max_value, NUMBER_OF_BINS+1) # +1 neede so number of bins is correct

    # calculate distribution with same binning
    hist_signal, _ = np.histogram(signal, bins)
    hist_background, _ = np.histogram(background, bins)

    # normalize histogram to 1 if needed
    if normalize:
        hist_signal = hist_signal/hist_signal.sum()
        hist_background = hist_background/hist_background.sum()
    
    return hist_signal, hist_background


def plot_separation(signal, background, variable):
    # create plot
    plt.figure()

    # needed weights to scale histogram to 1
    weights_signal = np.ones_like(signal)/float(len(signal))
    weights_background = np.ones_like(background)/float(len(background))

    # create histogram
    plt.hist(
        [signal, background],
        NUMBER_OF_BINS,
        weights=[weights_signal, weights_background], # normalize
        label=["Signal","Background"],
        histtype='step',
        stacked=False,
        fill=False
    )
    
    # beautificate
    plt.legend()
    
    min_value = np.min([np.min(signal), np.min(background)])
    max_value = np.max([np.max(signal), np.max(background)])
    plt.xlim([min_value,max_value])

    plt.savefig(f"separation_{variable.split('/')[-1]}.png")

    if SHOW_PLOTS:
        plt.show()
    else:
        plt.close()


def calculate_separation_power(signal, background):
    # create even histograms
    hist_signal, hist_background = transform_to_histogram(signal, background, True)

    # calculate separation power
    separation_power = 0

    for (bin_signal, bin_background) in zip(hist_signal, hist_background):
        # check for empty bins
        if(bin_signal+bin_background==0):
            continue
        
        separation_power+= (bin_signal-bin_background)**2 / (bin_signal+bin_background)

    # correctly scale separation power for double counted events
    separation_power/= 2

    return separation_power


if __name__=="__main__":
    main()