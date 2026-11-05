#!/usr/bin/env python3

"""Kinematic analysis of grasp and lift
Organization: The script requires the following organization of data:
    this script: put in directory with your subject folder data
    run data: both the *experiment_data file and all *motion/force tracking files in /Data/S**/ (example Data/S03/S03_experiment_data3217.txt)
    summary files will be saved in results folder: /results/S**/S**summary_starttrial*_.txt
    LoadPhase to Onset files will be saved in the results folder: /results/S**/S**trial**_lp_to_onset_starttrial*_.txt

    To run: open ipython, type "run Analysis_Script_Exp_Two3.py"

    While running, it waits for user input:
        1. At first plot, it waits for a 'y' confirmation, or accepts a mouse click to change the onset position, 
        or accepts left/right arrow to nudge onset one index backward/forward. Type 'y' to proceed

        2. After confirmation of first plot, it waits for a second 'y' confirmation, this time for roll. At this point,
        it also accepts left/right arrow to nudge roll one index backward/forward. 'y' to proceed
        
        3. After confirmation of the roll point, it waits for a third 'y' confirmation, this time for the valley after the max roll. 
        At this point it also accepts left/right arrow to nudge valley index backward/forward. 'y' to proceed

        4. Instead of 'y' confirmation, user can input 'b' to indicate a "bad" trial.  In this case, all relevant summary values 
        will be set to 0. And the script will proceed to the next trial

        5. After confirmation of the roll, user sees plot of the x,y,z positions and tangential velocities of thumb and pointer
        from the start of trial to a little after the lift onset.  It also shows the identified thumb and index reach onset, max tangential velocity,
        onset of thumb on T and onset of index on T and lift onset.
        Click 'c' to "continue" and move on.
        
        6. After the position and velocity plots, user sees grip force and lift force plots for the thumb and pointer finger. The force curves show 
        the whole data set for the trial. It also shows the lift force initiation and contact points for each digit. 
        Click 'c' to continue and move on. 
        
        7. After the digit force plots, user sees compensatory moment (Mcom) plots. The red dot is the point of lift onset. The plots ends around the endtone.
        Click 'c' to continue and move on.  

    This program calls numfxns.py,line_picker.py (custom) and detect_peaks.py (from https://github.com/demotu/BMC) - under /src directory

experiment_data file contents (columns): 

    0. Trial number -- good
    1. TR number at start -- good
    2. ITI -- good
    3. Center of Mass Condition (1 means left and 2 means right)
    4. Experimental condition (aim2: )
    5. Sync signal timestamp (TSYNC)
    6. First real trigger timestamp
    8. Timestamp of trigger at start of trial (self.triggered)
    9. Start tone timestamp
    10. End tone timestamp
    11. Time at which object is really down (self.treallydown)
    12. Start of intertrial interval
    13. End of intertrial interval
    14. Button-up timestamp
    15. Button-down timestamp
    
    
motion tracking file contents (columns):

    0. Timestamp
    1-3. X,Y,Z positions of the T
    4-6. Yaw, Pitch, Roll of the T
    7-9  X,Y,Z of thumb positions
    10-12 X,Y,Z of pointer positions


summary file contains (columns):

(0) totaltrialnum: trial number
(1) com_cond: center of mass condition (1 means left and 2 means right)
(2) firsttriggertime: first trigger time 
(3) starttone: start tone time 
(4) buttonuptime: time when participant releases button
(5) onsettime: lift onset timestamp
(6) buttondowntime: time when participant pushes the button after setting object down
(7) data['filter']['Tx'][onsetindx]: x position of the T at lift onset
(8) onsety: y position of the T at lift onset (equivalent to data['filter']['Ty'][onsetindx]) 
(9) data['filter']['Ty'][onsetindx]: y position of thr T at lift onset - note this has not taken into account what z was at start tone index like what is done for onsety
(10) data['filter']['Tz'][onsetindx]: z position of the T at lift onset
(11) data['filter']['rolloffset'][onsetindx]: roll at lift onset index
(12) absmaxroll: value of max roll within 250 ms window after lift onset
(13) onsetofmaxroll: time from lift onset to absmaxroll
(14) reachonset: time from first real trigger to button up time
(15) reachphase: time from button up time to first digit contact
(16) executeonset: time from first real trigger to lift onset time
(17) executephase: time from lift onset to button down time
(18) Contact_COP_difference: center of pressure difference at contact = thumb COP @ contact - pointer finger COP @ contact. Contact is defined when both thumb and pointer finger grip forces > 0.2 N. 
(19) gf_rate_thumb_max: maximum grip force rate of thumb from button up to 4 seconds later.
(20) gf_rate_pointer_max: maximum grip force rate of index finger from button up to 4 seconds later. 
(21) lf_rate_thumb_max: maximum lift force rate of thumb
(22) lf_rate_pointer_max: maximum lift force rate of index finger
(23) COP_diff_at_lift_onset: center of pressure difference at lift onset = thumb COP @ lift onset - pointer finger COP @ lift onset
(24) COP_pointer_at_lift_onset: pointer finger COP at lift onset
(25) COP_thumb_at_lift_onset: thumb COP at lift onset
(26) gf_thumb_at_lift_onset: thumb grip force at lift onset
(27) gf_pointer_at_lift_onset: pointer finger grip force at lift onset
(28) lf_thumb_at_lift_onset: thumb lift force at lift onset
(29) lf_pointer_at_lift_onset: pointer finger lift force at lift onset
(30) load_phase: Time from lift force initaition to lift onset
(31) Mcom_at_lift_onset: compensatory moment at lift onset
(32) Mcom[force_maxrollindx]: compensatory moment at max roll
(33) lf_diff_at_lift_onset: the lift force difference at lift onset 
(34) fbackminroll - min roll between 250ms after peak roll and the end tone 
(35) fbackmaxroll - max roll between 250ms after peak roll and the end tone 
(36) fbackmeanroll - mean roll between 250ms after peak roll and the end tone 
(37) fbackstdnroll - stdev roll between 250ms after peak roll and the end tone
(38) roll_at_valley: roll at valley after peak 
(39) roll_to_vall_time: time from from peak roll to first valley
(40) timeNet_lf_initiate
(41) timeNet_gf_initiate
(42) gfLFCoupling

Note: A roll to the left is a positive roll. A positive Mcom is a moment to the right.

Outputs: a summary txt file: (1) 

"""
__author__ = "Deborah Barany/Michelle Marneweck/Joshua Bland"
__version__ = "5.0 first implemented Oct 3, 2024"

import numpy as np
import sys
import matplotlib.pyplot as plt # type: ignore
import os
import glob
from detect_peaks import detect_peaks
from numfxns import fcount, butter_lowpass_filter,find_nearest_above, find_nearest
from line_picker import line_picker

#Constant Variable List
order = 4  #order of the butterworth filter (can try different orders)
fs_motion = 150.0 #sample rate in Hz
fs_force = 500.0 #sample rate in Hz
cutoff = 6 #desired cutoff frequency of the filter, Hz
downsample = 1  #if you wanted to skip some samples (to animate, for example), input a higher number
windowtime = .25  #250 ms
indxtime_motion = np.round(fs_motion*windowtime) #since refresh rate is 150 hz (150 samples/sec), we round up 
indxtime_force = np.round(fs_force*windowtime) #since refresh rate is 150 hz (150 samples/sec), we round up
sample_period_force = 1/fs_force
sample_period_motion = 1/fs_motion

#Name the tracking file columns
column_names_motion = {0:'mytime',1:'Tx',2:'Ty',3:'Tz',4:'Tyaw',5:'Tpitch',6:'Troll',7:'Thumb_x',8:'Thumb_y',9:'Thumb_z',10:'Pointer_x',11:'Pointer_y',12:'Pointer_z'}
force_column_names = {0:'time_stamps',1:'Fx32339',2:'Fy32339',3:'Fz32339',4:'Tx32339',5:'Ty32339',6:'Tz32339',7:'Fx32340',8:'Fy32340',9:'Fz32340',10:'Tx32340',11:'Ty32340',12:'Tz32340',13:'Arduino'}

#sets parent directory of your working directory. You should start this script in working directory, and there should be an aging/Data folders that follow this directory
DIR = '/Users/motorskill/Desktop/aging/'#DIR = os.path.normpath(os.getcwd() + os.sep + os.pardir) #'/Users/mirandahulsey-vincent'
print ('DIR',DIR)

#Input all subjects to analyze here (e.g., ["S01","S02","S03"])
subject_list = ["S1_blocked"]

# Set the trial_begin before running
#########################
trial_begin = 1    ######
#########################

# Variable letting the script know whether or not the starting COM was changed at least minute (after condition file created)
# Variable can be set to "yes" or "no"
#########################
COM_change =  "no" ######
#########################

# The grip width and grip surface thickness values are slightly different from experiment 1 to experiment 2 
grip_width = 46  
grip_surface_thickness = 3

summary = np.zeros((43,), dtype=int)

onsettime = 0

"""Function for rate calculations. Midpoint numerical differentiation is used to approximate instantaneous force rate at each data 
point. Midpoint method is used because it approximates the rate at each data point instead of halfway between data points."""
def rate_calculation_force(data):
	storage = []
	for i in range(0,len(data)):
		if i == 0:
				# One sided forward difference for first point
			instant_rate = (data[i+1] - data[i])/sample_period_force
			storage.append(instant_rate)
		elif i == (len(data)-1):
				# One sided backward difference for last point
			instant_rate = (data[i]-data[i-1])/sample_period_force
			storage.append(instant_rate)
		else:
				# Central difference for all other points
			instant_rate = (data[i+1] - data[i-1])/(2*sample_period_force)
			storage.append(instant_rate)
	return storage 

""" Function for velocity calculations. Midpoint numerical differnetiation used to approximate instantaneous velocity at each data 
point. Midpoint method is used because it approximates the velocity at each data point instead of halfway between data points. """
def velocity_calculation(x_data,y_data,z_data):
	storage = []
	for i in range(0,len(x_data)):
		if i == 0:
				# One sided forward difference for first point
			instant_rate = ((((x_data[i+1] - x_data[i])**2) + ((y_data[i+1] - y_data[i])**2) + ((z_data[i+1] - z_data[i])**2))**(0.5))/sample_period_motion
			storage.append(instant_rate)
		elif i == (len(x_data)-1):
				# One sided backward difference for last point
			instant_rate = ((((x_data[i] - x_data[i-1])**2) + ((y_data[i] - y_data[i-1])**2) + ((z_data[i] - z_data[i-1])**2))**(0.5))/sample_period_motion
			storage.append(instant_rate)
		else:
				# Central difference for all other points
			instant_rate = ((((x_data[i+1] - x_data[i-1])**2) + ((y_data[i+1] - y_data[i-1])**2) + ((z_data[i+1] - z_data[i-1])**2))**(0.5))/(2*sample_period_motion)
			storage.append(instant_rate)
	return storage 

#animation functions that need to be in main script
def press(event):
    print('press', event.key)
    sys.stdout.flush()
    N=event.key
    mutable_object['key'] = N

def onpick2(event):
    global mutable_object
    print(('onpick2 line:', event.pickx, event.picky))
    mutable_object['xdata'] = event.pickx
    mutable_object['ydata'] = event.picky
   
#Subject loop
for subject in subject_list:
    # each subject has own results folder. (Ex. /results/S01/)
    path=os.path.join(DIR,'data',subject)
    print('pathToSubjectData',path)
   
    #load experiment file name and path 
    exp_filename = glob.glob(os.path.join(path,subject+"*experiment_data"+"*"+"_"+str(trial_begin)+".txt"))
    print('expFileName',exp_filename)
  
    #exit out if script did not find an experiment file
    if len(exp_filename) != 1:
        assert 0

    exp_file = open(exp_filename[0],'r') # open the epxperiment file to read it.
    lines = exp_file.readlines() # every line of the file becomes an element in a list
    trials = len(lines) #number of trials
    starttrial = 1
    firstgoodtrial = 0 

    #loop through each trial in a given run's experiment file 
    for trialnum in range(starttrial,trials+1): 
        print('trialnum',trialnum)
        
        #Start with clear values
        absmaxroll = 0
        onsetofmaxroll = 0
        onsettime = 0
        onsety = 0
        firsttriggertime = 0
        buttonuptime = 0 
        buttondowntime = 0
        reachonset = 0
        reachphase = 0
        executeonset = 0
        executephase = 0
        Contact_COP_difference = 0
        gf_rate_thumb_max = 0
        gf_rate_pointer_max = 0
        lf_rate_thumb_max = 0
        lf_rate_pointer_max = 0
        Mcom_at_lift_onset = 0
        load_phase = 0
        COP_thumb_at_lift_onset = 0
        COP_pointer_at_lift_onset = 0
        COP_diff_at_lift_onset = 0
        lf_thumb_at_lift_onset = 0
        lf_pointer_at_lift_onset = 0
        gf_thumb_at_lift_onset = 0
        gf_pointer_at_lift_onset = 0
        Mcom_at_max_roll = 0 
        lf_diff_at_lift_onset = 0
        roll_at_valley = 0
        roll_to_vall_time = 0
        fbackmaxroll = 0
        fbackminroll = 0
        fbackmeanroll = 0
        fbackstdroll = 0
        gfLfCoupling = 0
        timeNet_lf_initiate = 0
        timeNet_gf_initiate = 0

        #actual trial number, given the trial the script was started at
        totaltrialnum = (trial_begin - 1) + trialnum    
        print('totaltrialnum', totaltrialnum)
        
        #The experiment data columns for the current run NEED TO CHANGE THESE BASED ON KEY TO COLUMNS
        exp_data = lines[trialnum-1].strip().split('\t')
        com_cond=float(exp_data[3]) #Center of Mass Condition (1 means left and 2 means right)
        print("com=",com_cond)
        tsyncsig = float(exp_data[5])   # Sync signal timestamp - time at which TTL voltage signal sent from Arduino to NI DAQ
        firsttrigger= float(exp_data[6])  #timestamp of first real trigger
        starttone = float(exp_data[8])  #the starttone
        endtone = float(exp_data[9])  #timestamp of endtone
        button_up = float(exp_data[13])  #timestamp for button release
        button_down = float(exp_data[14]) #timestamp for button press after object set down
        print("check button_down:", button_down)
        buttonuptime = button_up 
        buttondowntime = button_down
        firsttriggertime = firsttrigger
        reachonset = buttonuptime-firsttriggertime

        #get the motion tracking data file for that trial CHECK IF CHANGED NAMING CONVENTION
        track_file = glob.glob(os.path.join(path,subject+"_"+str(float(totaltrialnum))+"*motion_data"+"*"+str(trial_begin)+".txt"))
        track_data = np.loadtxt(track_file[0])
        print(track_file)

        #get the force transducer data file for that trial CHECK IF CHANGED NAMING CONVENTION
        force_file = glob.glob(os.path.join(path,subject+"_"+str(float(totaltrialnum))+"*force_data*"+str(trial_begin)+".txt"))
        print('force_file:',force_file)
        force_data = np.loadtxt(force_file[0], delimiter = ',')
        print(force_file)

        # The 13th column of the force data file are values from the Arduino input channel 
        arduino_signal = force_data[:,12]
        
        """Following timestamps for force transducer data using voltage signal from Arduino to DAQ. 
        The time at which TTL voltage signal (tsync signal) is sent from the Arduino to the DAQ is known - 
        this is matched with the data point at which the Arduino voltage input channel on the DAQ jumps above 5 volts. 
        The other timestamps of the force transducer data are worked out using the known sampling rate. Maximum time
        it could be off is less than 0.002 seconds based on sampling rate."""
        
        for volt_point in range(0,len(arduino_signal)): #Find where the Arduino voltage inout channel is first > 5v
            if arduino_signal[volt_point] > 4:
                sync_point = volt_point 
                break  
        print("sync point:",sync_point)
        initial_force_time = tsyncsig - (sync_point*0.002) #Find timestamp of first force transducer data point
        timestamps = np.linspace(initial_force_time,(initial_force_time+((len(arduino_signal)-1)*0.002)),4000) #Array of timestmaps for all force transducer data points
        print("initial_force_time:",initial_force_time)
        timestamps2 = timestamps[np.newaxis].T   
        force_data_2 = np.hstack((timestamps2,force_data)) # Append timestamp column to the rest of forcr transducer data
        print(force_data_2.shape)
        
        
        #create dictionary of the force data           
        data_force = {}
        data_force['columns'] = force_column_names 
        data_force['raw'] = {}
        data_force['filter'] = {}
    
        # Filter the force data
        for col in range(len(data_force['columns'])):
            raw_data = force_data_2[:,col][::downsample]
            data_force['raw'][data_force['columns'][col]] = raw_data
            data_force['filter'][data_force['columns'][col]] = butter_lowpass_filter(raw_data,cutoff,fs_force,order)
        time_stamps = data_force['raw']['time_stamps']

        """Following creates arrays for grip force and lift force at every data point. Every single data point is subtracted
        from the average of the first 250 data points to account for drift over trials - the force sensors are not being touched 
        during the first 250 data points and so this average acts as a bias offset. The absolute value is taken for the normal force 
        in order to reorient the coordinate system so that grip forces are positive."""
        gf_thumb = abs(data_force['filter']['Fz32339']) - np.mean(abs(data_force['filter']['Fz32339'][:250])) # Thumb grip force
        gf_pointer = abs(data_force['filter']['Fz32340']) - np.mean(abs(data_force['filter']['Fz32340'][:250]))  # Pointer finger grip force
        lf_thumb = data_force['filter']['Fy32339'] - np.mean(data_force['filter']['Fy32339'][:250])  # Thumb lift force
        lf_pointer = data_force['filter']['Fy32340'] - np.mean(data_force['filter']['Fy32340'][:250]) # Pointer finger lift force
            
        lf_diff = lf_thumb - lf_pointer #Lift force difference (thumb - pointer finger)
        print('Tx32339 units check:', np.mean(data_force['filter']['Tx32339'][:100]))
        #Following calcualtes grip force rates using the "rate_calculation_force" equation defined at the beginning of script. 
        gf_rate_thumb = rate_calculation_force(gf_thumb) #grip force rate at every data point for thumb
        gf_rate_pointer = rate_calculation_force(gf_pointer) # grip force rate for every data point for pointer finger
        
        #Following calculates lift force rates using the "velocity_calculation" equation defined at the beginning of script. 
        lf_rate_thumb = rate_calculation_force(lf_thumb) # lift force rate for every data point for for thumb
        lf_rate_pointer = rate_calculation_force(lf_pointer) #lift force rate for every data point for pointer finger
        
        
        # Finding the force data timestamp index for nearest neighbor of the endtone 
        set_down_force_indx = find_nearest(time_stamps,endtone)

        #find force data timstamp index for nearest neighbor of startTone index
        startToneForceindx=find_nearest(time_stamps,starttone)
        print('startToneForceindx',startToneForceindx)
        
        #find the force data timestap index for nearest neighbor of buttonup 
        buttonupForceindx=find_nearest(time_stamps,button_up)
        print('buttonupForceindx',buttonupForceindx)
        
        # Maximum grip force and lift force rates from start of trial to around when the object should be set down
        gf_rate_thumb_max = max(gf_rate_thumb[0:(set_down_force_indx+1)]) 
        gf_rate_pointer_max = max(gf_rate_pointer[0:(set_down_force_indx+1)])
        lf_rate_thumb_max = max(lf_rate_thumb[0:(set_down_force_indx+1)])
        lf_rate_pointer_max = max(lf_rate_pointer[0:(set_down_force_indx+1)])

        # Finding point of contact for thumb. Defined as point when exceeds 0.2 N threshold. 
        contact_point_thumb = next((i for i in range(buttonupForceindx, len(gf_thumb)) if gf_thumb[i] > 0.2), -999)
        print('contact_point_T', contact_point_thumb)
        contact_point_pointer = next((i for i in range(buttonupForceindx, len(gf_pointer)) if gf_pointer[i] > 0.2), -999)
        print('contact_point_P', contact_point_pointer)

        if contact_point_thumb == -999:
        # Finding point of contact for thumb. Defined as point when exceeds 0.2 N threshold. 
            contact_point_thumb = next((i for i in range(startToneForceindx, len(gf_thumb)) if gf_thumb[i] > 0.2), -999)
            print('contact_point_T from startToneindx found', contact_point_thumb)
        
        if contact_point_pointer == -999:
            contact_point_pointer = next((i for i in range(startToneForceindx, len(gf_pointer)) if gf_pointer[i] > 0.2), -999)
            print('contact_point_P from startToneindx found', contact_point_pointer)

        gf_contact = -999
        reachphase = -999

        if contact_point_thumb != -999 or contact_point_pointer != -999:
            if contact_point_thumb >= buttonupForceindx and contact_point_pointer >= buttonupForceindx:
                gf_contact = min(contact_point_thumb, contact_point_pointer)
                print('gfcontact',gf_contact)
                first_contact = gf_contact
                print('firstContact', first_contact)
                reachphase = time_stamps[first_contact]-buttonuptime
                print('reachphase',reachphase)
            if contact_point_thumb >= startToneForceindx and contact_point_pointer >= startToneForceindx:
                gf_contact = min(contact_point_thumb, contact_point_pointer)
                print('gfcontact calculated from startTone',gf_contact)
                first_contact = gf_contact
                print('firstContact calculated from startTone', first_contact)
                reachphase = time_stamps[first_contact]-buttonuptime

        """x-axis torque values for the thumb and pointer finger. Every single data point is subtracted
            from the average of the first 250 data points to account for drift over trials - the force sensors 
            are not being touched during the first 250 data points and so this average acts as a bias offset."""           
        tx_thumb = (data_force['filter']['Tx32339'] - np.mean(data_force['filter']['Tx32339'][:250]))*1000
        tx_pointer = (data_force['filter']['Tx32340'] - np.mean(data_force['filter']['Tx32340'][:250]))*1000

        """ Center of pressure values for thumb and pointer finger at every data point. Center of pressure is vertical 
        distance from the center of rhte force sensor. """
        COP_thumb = (-tx_thumb - lf_thumb*grip_surface_thickness)/gf_thumb
        COP_pointer = ((-tx_pointer) - (lf_pointer*grip_surface_thickness))/gf_pointer
        COPdifference = COP_thumb - COP_pointer # Center of pressure difference at evey data point (thumb - pointer)
        print('Thimb tX', tx_thumb)
        Contact_COP_difference = COP_thumb[gf_contact] - COP_pointer[gf_contact] # Center of pressure difference at 'overall' contact point.

        if gf_contact == -999:
              Contact_COP_difference = -999
            
        #Determining load force initiation for thumb. Defined as point when exceeds 0.2 N threshold and keeps increasing for 20 samples.  
        net_lift_force = lf_thumb + lf_pointer 
        mean_grip_force=(gf_thumb+gf_pointer)/2

        """ The following determines the 'overall' load force initiation index.  The 'overall' load force initiation 
        point is defined as when both the thumb and the pointer finger have reached their load force initiation points. """
            
        for i in range(0,len(net_lift_force)): #goes through the thumb lift force array.
            if net_lift_force[i] > 0.2 and np.all(np.diff(net_lift_force[i:i+21]) > 0): # once the threshold is passed and it increseas after, record index & break loop 
                net_lf_initiate = i
                found_net_lf_initiate = True
                break
            else:
                found_net_lf_initiate = False
                net_lf_initiate = len(net_lift_force)-1
                
        for i in range(0,len(mean_grip_force)): #step to find mean grip force initiation - goes through the thumb grip force array.
            if mean_grip_force[i] > 0.2 and np.all(np.diff(mean_grip_force[i:i+20]) > 0): # once the threshold is passed and it increseas after, record index & break loop 
                mean_gf_initiate = i
                found_mean_gf_initiate = True
                break
            else:
                found_mean_gf_initiate = False
                mean_gf_initiate = len(mean_grip_force)-1
                
        gfLfCoupling = time_stamps[net_lf_initiate] - time_stamps[mean_gf_initiate]
        timeNet_lf_initiate= time_stamps[net_lf_initiate]
        timeNet_gf_initiate=time_stamps[mean_gf_initiate]
        
        #create dictionary of the motion tracking data            
        data = {}
        data['columns'] = column_names_motion
        data['raw'] = {}
        data['filter'] = {}
        
        #extract and filter (smooth out) motion tracking data             
        for col in range(len(data['columns'])):
            raw_data = track_data[:,col][::downsample]
            data['raw'][data['columns'][col]] = raw_data
            # if col == 6: #because roll is column 7   
            #     raw_data = np.abs(raw_data)
            data['filter'][data['columns'][col]] = butter_lowpass_filter(raw_data,cutoff,fs_motion,order)
        mytime = data['raw']['mytime']

        #data['raw']['rolloffset']=(np.abs(data['raw']['Troll']))
        data['filter']['rolloffset'] = data['filter']['Troll']

        # # Assuming Tx_thumb and gf_thumb are arrays or lists of equal length
        # for i in range(len(tx_thumb)):
        #     if gf_thumb[i] != 0:
        #         cop_val = abs(tx_thumb[i] / gf_thumb[i])
        #         print(f'Trial {i}: abs(Tx_thumb / gf_thumb) = {cop_val:.2f}')
        #     else:
        #         print(f'Trial {i}: gf_thumb is zero — cannot compute COP')
        #index of the start tone for motion tracking data
        startindx = find_nearest(mytime,starttone) 
        endindx=find_nearest(mytime,endtone)

        # Switching the center of mass condition if COM_change == yes
        if com_cond == 1 and COM_change == "yes":   
            com_cond = 2
        elif com_cond == 2 and COM_change == "yes":
            com_cond = 1      
        
        ##subtract or add roll and object Y values at start tone from roll and y position of object values  
        if data['filter']['rolloffset'][startindx]>0:
            data['filter']['rolloffset']=data['filter']['rolloffset']-data['filter']['rolloffset'][startindx]
        else:
            data['filter']['rolloffset']=data['filter']['rolloffset']-data['filter']['rolloffset'][startindx]            
        if data['filter']['Ty'][startindx] <0:
           data['filter']['objectyoffset']=data['filter']['Ty']+(np.abs(data['filter']['Ty'])[startindx])
        else:
            data['filter']['objectyoffset']=data['filter']['Ty']-(np.abs(data['filter']['Ty'])[startindx])
            
        ###FIND LIFT ONSET
        #Find all valleys and peaks using detect_peaks function 
        valleysy = detect_peaks(data['filter']['objectyoffset'],mpd=1,edge='rising',valley=True)
        peaksy = detect_peaks(data['filter']['objectyoffset'],mpd=1,edge='rising',valley=False)

        #Find the onset based on peaks and valleys
        foundonset = 0
        onsetindx = 0
        trialsummary=1
        isBadOnset = False
        
        #if the max y value of the object is smaller or equal to the y value of the object at the start index then it's an indication that object wasn't lifted, so flag this as bad trial
        if np.max(data['filter']['objectyoffset'][startindx:]) <=  data['filter']['objectyoffset'][startindx]:
            trialsummary=0
            print('tracking bad 1')
            isBadOnset = True
            
        #flag trial as no good if max y value of object is equal to 1cm or less from start index and beyond
        if np.max(data['filter']['objectyoffset'][startindx:]) <=  0.005:
            trialsummary=0
            print('tracking bad 2')
            isBadOnset = True
         
        if isBadOnset == False:
        #tries to find the lift onset - opens up a visualization of the lift and jumps to point where it thinks lift onset happened. Allow user to toggle to correct/check lift onset
        #find valley where the next peak is much* (.1) bigger than the valley
            for valley in valleysy: #cycle through valleys
                nextpeak = find_nearest_above(peaksy,valley)
                #add constraint that peak can't be right next to valley:
                if peaksy[nextpeak] - valley:
                    nextpeak = nextpeak + 1
                    print('nextpeak',nextpeak)
                #check if next peak is too high, and flag it
                if len(peaksy) == nextpeak:
                    trialsummary=0
                    print('next peak too high')
                    break
                #find the valley in which there's a change aka lift by the next peak
                if data['filter']['objectyoffset'][peaksy[nextpeak]] > data['filter']['objectyoffset'][valley] + .005:  #while next peak is 5mm bigger than the valley before it, set onset to on and onsetindx to valley
                    foundonset = 1
                    onsetindx = valley
                    #refine estimate of valley so that y is above 1 mm
                    #loop until found valley > 0.001 and make sure deriative of the 20 values from onsetindex are all above 0 (i.e., values are increasing)
                    print('TOTAL LEN: ',len(data['filter']['objectyoffset']))
                    while len(data['filter']['objectyoffset']) > onsetindx + 20 and (data['filter']['objectyoffset'][onsetindx] < 0.001 or np.all(np.diff(data['filter']['objectyoffset'][onsetindx:onsetindx+20])>=0)==False): #for as long valueatonsetindx is less than 1mm or there arent 20 samples bigger than 0
                        onsetindx = onsetindx + 1
                    else:
                        print('BREAK on onsetindx',onsetindx)
                        break

        #few more final checks of bad onsets and will then give up if no onset was found.        
        if trialsummary == 0:
            print('continuing because could not find good onset')
            isBadOnset = True
           
        if np.sum(np.diff(data['filter']['objectyoffset'][onsetindx:onsetindx+20])>=0)/float(len(data['filter']['objectyoffset'][onsetindx:onsetindx+20])) < .8:
            trialsummary=0
            print('tracking bad 3')
            isBadOnset = True
        
        if onsetindx+1 == len(data['filter']['objectyoffset']):
            print('bad onset finding')
            trialsummary = 0
            isBadOnset = True
        if foundonset == 0:  #didn't find onset
            print('DID NOT FIND ONSET, SETTING TO 0')
            onsetindx = 0
            isBadOnset = True
        
        #get initial onset information
        onset = mytime[onsetindx] #onset time
        start=mytime[startindx]
        print('mytimestart',start)
        print('mytimeatonset',onset) 
        onsetybeforesubtraction=data['filter']['Ty'][onsetindx]
        onsety = data['filter']['objectyoffset'][onsetindx]  #y position of tracker at onset
        yatstart = data['filter']['objectyoffset'][startindx]
        yatstartbeforesubtraction=data['filter']['Ty'][startindx]

        #get the roll window from lift onset to 250ms later
        rollwindow = data['filter']['rolloffset'][onsetindx:onsetindx+int(indxtime_motion)]

        #find initial absolute value and timepoint of max roll within 250 ms window after onset
        absmaxroll = np.max(np.abs(rollwindow))
        absmaxrollindx = np.argmax(np.abs(rollwindow))
        print('absolute value of max roll', absmaxroll)
        
        #time relative to onset when max roll was found
        onsetofmaxroll = mytime[onsetindx+absmaxrollindx] - mytime[onsetindx] 
        mytimemaxroll = mytime[onsetindx+absmaxrollindx]
        print('mytimemaxroll', mytimemaxroll)
        found_valley = True
        
        if (data['filter']['rolloffset'][(onsetindx+absmaxrollindx)])>0 and ((onsetindx + absmaxrollindx) < endindx):
            valleysr = detect_peaks(data['filter']['rolloffset'][(onsetindx+absmaxrollindx):endindx],mpd=1,edge='rising',valley=True)
            if len(valleysr) == 0:
                found_valley = False 
            else:
                valley_index = onsetindx+absmaxrollindx+valleysr[0]
        elif (data['filter']['rolloffset'][(onsetindx+absmaxrollindx)]) < 0 and ((onsetindx + absmaxrollindx) < endindx):
            valleysr = detect_peaks(data['filter']['rolloffset'][(onsetindx+absmaxrollindx):endindx],mpd=1,edge='rising',valley=False)
            if len(valleysr) == 0:
                found_valley = False 
            else:
                valley_index = onsetindx+absmaxrollindx+valleysr[0]    
        else: 
            found_valley = False
            
        #Visually confirm the correct start, roll and valley after peak roll has been found 
        for check in ['onset','roll','valley']:
            confirm = 'n'
            count = 0
            if check == 'onset':
                onsetchange = 0
                print('check onset')
                print('com_cond',com_cond)
            elif check == 'roll':
                rollchange = 0
                print('check roll')
            elif check == 'valley':
                valleychange = 0
                print('check valley')      
            mutable_object = {} 

            fig = None
            while confirm != 'y' and confirm != 'b' and isBadOnset == False:
                #plot initial onset and roll information
                if count == 0:
                    fig, ax = plt.subplots(2)
                    count = 1

                    #accepts user input about key press
                    fig.canvas.mpl_connect('key_press_event', press)

                    #plot y position of T on upper subplot and accept user input for mouse click on line
                    line, = ax[0].plot(mytime,data['filter']['objectyoffset'], picker = line_picker)
                    fig.canvas.mpl_connect('pick_event', onpick2)
                    fig.tight_layout(h_pad=1.5)
                    plt.subplots_adjust(top=0.89)
                    if onsetchange == 0 or check == 'roll':
                        point = ax[0].scatter(onset,onsety,c='r')
                    ax[0].set_title('Y Position of T')

                    #plot roll information on bottom subplot
                    ax[1].plot(mytime,(data['filter']['rolloffset']))
                    fig.suptitle('trial number '+str(trialnum),fontweight="bold")
                    line1 = ax[1].scatter(np.ones(20)*onset,np.linspace(np.min(data['filter']['rolloffset']),np.max(data['filter']['rolloffset']),20),marker = '|')
                    line2 = ax[1].scatter(np.ones(20)*(onset+windowtime),np.linspace(np.min(data['filter']['rolloffset']),np.max(data['filter']['rolloffset']),20),marker = '|')
                    ax[1].set_title('Roll of T')
                    
                    if check == 'onset' or check == 'valley' or (check == 'roll' and rollchange == 0):  
                        rollpoint = ax[1].scatter(mytime[onsetindx+absmaxrollindx],(data['filter']['rolloffset'][onsetindx+absmaxrollindx]))
                    if check == 'onset' or check == 'roll' or (check == 'valley' and valleychange == 0):
                        if found_valley:
                            valleypoint = ax[1].scatter(mytime[valley_index],(data['filter']['rolloffset'][valley_index]),c='b')
                        #rollpoint_a = ax[1].scatter(mytime[onsetindx+absmaxrollindx_a],data['filter']['Troll'][onsetindx+absmaxrollindx_a],c='y')             
                    plt.draw()
                plt.pause(.2) 

                #If user clicked or pressed key
                if 'key' in mutable_object or 'xdata' in mutable_object:
                    #get key press information, if any
                    if 'key' in mutable_object:
                        confirm = mutable_object['key']
                    #get any information that user wants to change the onset or roll values, and update accordingly
                    if 'xdata' in mutable_object or mutable_object['key'] == 'right' or mutable_object['key'] == 'left':
                        line1.remove()
                        line2.remove()
                        rollpoint.remove()
                        point.remove()
                        if found_valley:
                            valleypoint.remove()
                        if check == 'onset':
                            onsetchange = 1
                        elif check == 'roll':
                            rollchange = 1
                        elif check == 'valley':
                            valleychange = 1
                            
                        if 'key' in mutable_object:
                            if mutable_object['key'] == 'right':  #move marker one index forward
                                if check == 'onset':
                                    onsetindx = onsetindx + 1
                                elif check == 'roll':
                                     absmaxrollindx = absmaxrollindx + 1
                                     #absmaxrollindx_a = absmaxrollindx_a + 1
                                     absmaxroll = data['filter']['rolloffset'][onsetindx+absmaxrollindx]
                                     print('max roll', absmaxroll) 
                                elif check == 'valley' and found_valley:
                                    if valley_index == endindx-1:
                                        print('max increase for valley point - too close to button down time')
                                    else:
                                        valley_index = valley_index + 1
                                    
                            elif mutable_object['key'] == 'left': #move marker one index backward
                                if check == 'onset':
                                    onsetindx = onsetindx - 1
                                elif check == 'roll':
                                    absmaxrollindx = absmaxrollindx - 1
                                    #absmaxrollindx_a = absmaxrollindx_a - 1
                                    absmaxroll = data['filter']['rolloffset'][onsetindx+absmaxrollindx]
                                    print('max roll', absmaxroll)
                                elif check == 'valley' and found_valley:
                                    if valley_index == (onsetindx+absmaxrollindx)+1:
                                        print('max decrease for valley point - too close to max roll point')
                                    else:
                                        valley_index = valley_index - 1
                            
                            onset = mytime[onsetindx] #onset time
                            onsety = data['filter']['objectyoffset'][onsetindx]  #y position of tracker at onset
                            point = ax[0].scatter(onset,onsety,c='g')
                            rollpoint = ax[1].scatter(mytime[onsetindx+absmaxrollindx],(data['filter']['rolloffset'][onsetindx+absmaxrollindx]))
                            #rollpoint = ax[1].scatter(mytime[onsetindx+absmaxrollindx_a],data['filter']['Troll'][onsetindx+absmaxrollindx_a])
                            mutable_object = {}
                        else:
                            onset = mutable_object['xdata'][0][0]
                            onsety = mutable_object['ydata'][0][0]
                            point = ax[0].scatter(onset,onsety,c='b')
                            onsetindx = np.where(mytime == mutable_object['xdata'][0][0])[0][0]
                        line1 = ax[1].scatter(np.ones(20)*onset,np.linspace(np.min(data['filter']['rolloffset']),np.max(data['filter']['rolloffset']),20),marker = '|')
                        line2 = ax[1].scatter(np.ones(20)*(onset+windowtime),np.linspace(np.min(data['filter']['rolloffset']),np.max(data['filter']['rolloffset']),20),marker = '|')

                        if check == 'onset':
                             rollpoint = ax[1].scatter(mytime[onsetindx+absmaxrollindx],(data['filter']['rolloffset'][onsetindx+absmaxrollindx]))  
                           # rollpoint = ax[1].scatter(mytime[onsetindx+absmaxrollindx_a],data['filter']['Troll'][onsetindx+absmaxrollindx_a])
                             print('onset y value', onsety)
                        if check == 'onset' or check =='roll':
                            if (data['filter']['rolloffset'][(onsetindx+absmaxrollindx)])>0 and ((onsetindx + absmaxrollindx) < endindx):
                                valleysr = detect_peaks(data['filter']['rolloffset'][(onsetindx+absmaxrollindx):endindx],mpd=1,edge='rising',valley=True)
                                if len(valleysr) == 0:
                                    found_valley = False 
                                else:
                                    found_valley = True
                                    valley_index = onsetindx+absmaxrollindx+valleysr[0]
                            elif (data['filter']['rolloffset'][(onsetindx+absmaxrollindx)])< 0 and ((onsetindx + absmaxrollindx) < endindx):
                                valleysr = detect_peaks(data['filter']['rolloffset'][(onsetindx+absmaxrollindx):endindx],mpd=1,edge='rising',valley=False)
                                if len(valleysr) == 0:
                                    found_valley = False 
                                else:
                                    found_valley = True
                                    valley_index = onsetindx+absmaxrollindx+valleysr[0]
                            else:
                                found_valley = False
                                
                            if found_valley:
                                valleypoint = ax[1].scatter(mytime[valley_index],(data['filter']['rolloffset'][valley_index]),c='b')
                        if check == 'valley' and found_valley:
                            valleypoint = ax[1].scatter(mytime[valley_index],(data['filter']['rolloffset'][valley_index]),c='b')
                            new_valley_val = data['filter']['rolloffset'][valley_index]
                            print('valley point',new_valley_val)

                        if mytime[onsetindx] != onset or onsety != data['filter']['objectyoffset'][onsetindx]:
                            assert 0
                else:
                    confirm = 'n' #change back to n

                #reject mouse click if selected onset y position was below 0
                if onsety < 0:
                    point.remove()
                    point = ax[0].scatter(onset,onsety, c='m')
                    confirm = 'n'
                    print('select again')      

                #update roll information
                rollwindow = data['filter']['rolloffset'][onsetindx:onsetindx+int(indxtime_motion)]
                absmaxroll = data['filter']['rolloffset'][onsetindx+absmaxrollindx]
                onsetofmaxroll = mytime[onsetindx+absmaxrollindx] - mytime[onsetindx] 

                if not found_valley:
                    fbackmaxroll = -999
                    fbackminroll = -999
                    fbackmeanroll = -999
                    fbackstdroll = -999
                    roll_at_valley = -999
                    roll_to_vall_time = -999           
                else: 
                    fbackwindow= data['filter']['rolloffset'][valley_index:endindx]
                    fbackmaxroll = np.max(np.abs(fbackwindow))
                    fbackminroll = np.min(np.abs(fbackwindow))
                    fbackmeanroll = np.mean(np.abs(fbackwindow))
                    fbackstdroll = np.std(np.abs(fbackwindow))
                    roll_at_valley = data['filter']['rolloffset'][valley_index]
                    roll_to_vall_time = mytime[valley_index] - mytime[onsetindx+absmaxrollindx]
            if fig is not None:
                plt.close(fig)

        if confirm != 'b' and isBadOnset == False:          
            """Finding the force data timestamp index for nearest neighbor of the lift onset timestamp (from motion data). 
            Nearest neighbor could be a maximum of 0.004 seconds away due to difference in sampling frequency."""
            force_onsetindx = find_nearest(time_stamps,onset)
            force_maxrollindx = find_nearest(time_stamps,mytime[onsetindx+absmaxrollindx])
            
            # calculating the compensatory moment for every data point. 
            Mcom = (lf_diff * (grip_width/2)) + (((gf_thumb+gf_pointer)/2)*COPdifference)

            # calculate key force parameters at lift onset defined by force_onsetindx
            Mcom_at_lift_onset = Mcom[force_onsetindx]
            print("mcom at lift onset:", Mcom_at_lift_onset)
            Mcom_at_max_roll = Mcom[force_maxrollindx]
            COP_thumb_at_lift_onset = COP_thumb[force_onsetindx] 
            COP_pointer_at_lift_onset = COP_pointer[force_onsetindx]
            COP_diff_at_lift_onset = COP_thumb[force_onsetindx] - COP_pointer[force_onsetindx]
            lf_thumb_at_lift_onset = lf_thumb[force_onsetindx]
            lf_pointer_at_lift_onset = lf_pointer[force_onsetindx]
            gf_thumb_at_lift_onset = gf_thumb[force_onsetindx]
            gf_pointer_at_lift_onset = gf_pointer[force_onsetindx]
            lf_diff_at_lift_onset = lf_diff[force_onsetindx]
                     
            onsettime = onset
            executeonset = onsettime-firsttriggertime
            executephase = buttondowntime-onsettime

            if not found_net_lf_initiate:       
                load_phase = -999
                lplf_slice_thumb = -999
                lplf_slice_pointer = -999
                lpgf_slice_thumb = -999
                lpgf_slice_pointer = -999
                lplf_rate_slice_thumb = -999
                lplf_rate_slice_pointer = -999
                lpgf_rate_slice_thumb = -999
                lpgf_rate_slice_pointer = -999
                gfLfCoupling = -999
                timeNet_lf_initiate= -999
                timeNet_gf_initiate=-999  
            else:
                load_phase = onset - time_stamps[net_lf_initiate] #load phase is time from "overall" load force initiate to lift onset
                lplf_slice_thumb = lf_thumb[net_lf_initiate:force_onsetindx] 
                lplf_slice_pointer = lf_pointer[net_lf_initiate:force_onsetindx]
                lpgf_slice_thumb = gf_thumb[net_lf_initiate:force_onsetindx]
                lpgf_slice_pointer = gf_pointer[net_lf_initiate:force_onsetindx] 
                lplf_rate_slice_thumb =  lf_rate_thumb[net_lf_initiate:force_onsetindx]
                lplf_rate_slice_pointer = lf_rate_pointer[net_lf_initiate:force_onsetindx]
                lpgf_rate_slice_thumb = gf_rate_thumb[net_lf_initiate:force_onsetindx]
                lpgf_rate_slice_pointer = gf_rate_pointer[net_lf_initiate:force_onsetindx]          
            
            # Plot to verify contact point from grip force values and lift force initiation point 
            fig, ax = plt.subplots(2)
            fig.canvas.mpl_connect('key_press_event',press)
            while mutable_object['key'] != 'c':
                ms = 20 # marker size
                 # Subplot 1: lift force - thumb and pointer finger
                ax[0].set_title('Lift Force - Thumb (black) & Index (green)')
                line_thumblf = ax[0].plot(time_stamps,lf_thumb, c='k')
                line_pointerlf = ax[0].plot(time_stamps, lf_pointer, c='g')
                line_netlf = ax[0].plot(time_stamps,net_lift_force, c='r')
                line_netlf_initiate = ax[0].scatter(time_stamps[net_lf_initiate],net_lift_force[net_lf_initiate], c='b')
                 # Subplot 2: grip force - thumb and pointer finger
                ax[1].set_title('Grip Force - Thumb (black) & Index (green)')
                line_thumbgf = ax[1].plot(time_stamps,gf_thumb, c='k')
                line_pointergf =  ax[1].plot(time_stamps,gf_pointer, c='g')
                point_thumbgf_contact = ax[1].scatter(time_stamps[contact_point_thumb],gf_thumb[contact_point_thumb], c='k',s=ms)
                point_pointergf_contact = ax[1].scatter(time_stamps[contact_point_pointer],gf_pointer[contact_point_pointer], c='r',s=ms)
                 #PLOT GRIP FORCE ONCTACT POINTE HERE 
                plt.suptitle('trial number '+str(trialnum),fontweight="bold")
                plt.tight_layout()
                plt.draw()
                plt.pause(.5)
            plt.close(fig)                                  
            
            mutable_object['key'] = 'i'
            
            fig = plt.figure()
            fig.canvas.mpl_connect('key_press_event',press)
            while mutable_object['key'] != 'c':
                #plt.plot(time_stamps[gf_contact:set_down_force_indx],Mcom[gf_contact:set_down_force_indx], c ='k')
                plt.plot(time_stamps[buttonupForceindx:set_down_force_indx],Mcom[buttonupForceindx:set_down_force_indx], c ='k')
                plt.scatter(time_stamps[force_maxrollindx],Mcom[force_maxrollindx], c = 'b')
                plt.scatter(time_stamps[force_onsetindx],Mcom[force_onsetindx], c = 'r')
                if com_cond == 1.0:
                    plt.axhline(y=300, color='r', linestyle='--')
                elif com_cond == 2.0:
                    plt.axhline(y=-300, color='r', linestyle='--')                        
                plt.suptitle('trial number '+str(trialnum),fontweight="bold")
                plt.title('Mcom')
                plt.draw()
                plt.pause(.5)
            plt.close(fig)
            
            mutable_object['key'] = 'i'
            
            fig, ax = plt.subplots(3, sharex = True)
            fig.canvas.mpl_connect('key_press_event',press)
            while mutable_object['key'] != 'c':
                ax[0].axhline(y=0, color='k', linestyle='-', linewidth = 0.5)
                ax[0].plot(mytime,data['filter']['objectyoffset'], c ='k')
                ax[0].set_xlim([starttone,endtone])
                ax[0].set_title('Object Height')
                ax[0].axvline(x=time_stamps[first_contact], linestyle='dotted', c='b', linewidth =2)
                ax[0].axvline(x=onset, linestyle='dotted', c='g', linewidth =2)
                ax[1].axhline(y=0, color='k', linestyle='-', linewidth = 0.5)
                ax[1].plot(mytime,data['filter']['rolloffset'], c = 'k')
                ax[1].set_title('Object Roll')
                ax[1].axvline(x=onset, linestyle='dotted', c='g', linewidth =2)
                ax[1].axvline(x=time_stamps[first_contact], linestyle='dotted', c='b', linewidth =2)
                ax[1].axvline(x=buttonuptime, linestyle='dotted', c='k', linewidth =2)
                ax[2].axhline(y=0, color='k', linestyle='-', linewidth = 0.5)
                ax[2].plot(time_stamps[first_contact:set_down_force_indx],Mcom[first_contact:set_down_force_indx], c = 'k')
                if com_cond == 1.0:
                    ax[2].axhline(y=236.29, color='r', linestyle='--')
                elif com_cond ==2.0:
                    ax[2].axhline(y=-236.29, color='r', linestyle='--')                        
                ax[2].set_title('Mcom')
                ax[2].axvline(x=onset, linestyle='dotted', c='g', linewidth =2)
                ax[2].axvline(x=time_stamps[first_contact], linestyle='dotted', c='b', linewidth =2)
                fig.suptitle('trial number '+str(trialnum),fontweight="bold")
                plt.tight_layout()
                plt.draw()
                plt.pause(.5)
            plt.close(fig)      

          # Checking if the results folder exists. If not create one.      
        if os.path.exists(os.path.join(DIR,'aging/results',subject)) is False:
            os.makedirs(os.path.join(DIR,'aging/results',subject))
            print('Created New Folder')
    
        if isBadOnset == True or confirm == 'b':
            onsettime = -999
            data['filter']['Tx'][onsetindx]  = -999
            onsety = -999
            data['filter']['Ty'][onsetindx] = -999
            data['filter']['Tz'][onsetindx] = -999
            data['filter']['rolloffset'][onsetindx] = -999
            absmaxroll = -999
            onsetofmaxroll = -999
            executeonset = -999
            executephase = -999
            COP_diff_at_lift_onset = -999
            COP_pointer_at_lift_onset = -999
            COP_thumb_at_lift_onset = -999
            gf_thumb_at_lift_onset = -999
            gf_pointer_at_lift_onset = -999
            lf_thumb_at_lift_onset = -999
            lf_pointer_at_lift_onset = -999
            load_phase = -999
            Mcom_at_lift_onset = -999
            Mcom_at_max_roll = -999
            lf_diff_at_lift_onset = -999
            fbackminroll = -999
            fbackmaxroll = -999
            fbackmeanroll = -999
            fbackstdroll = -999
            roll_at_valley = -999
            roll_to_vall_time = -999
        print('absmaxroll',absmaxroll)
        trialsummary = np.array([totaltrialnum, com_cond, firsttriggertime, starttone, buttonuptime, onsettime, buttondowntime, data['filter']['Tx'][onsetindx], onsety, data['filter']['Ty'][onsetindx], data['filter']['Tz'][onsetindx], data['filter']['rolloffset'][onsetindx], absmaxroll, onsetofmaxroll, reachonset, reachphase, executeonset, executephase, Contact_COP_difference, gf_rate_thumb_max, gf_rate_pointer_max, lf_rate_thumb_max, lf_rate_pointer_max, COP_diff_at_lift_onset, COP_pointer_at_lift_onset, COP_thumb_at_lift_onset, gf_thumb_at_lift_onset, gf_pointer_at_lift_onset, lf_thumb_at_lift_onset, lf_pointer_at_lift_onset,load_phase, Mcom_at_lift_onset,Mcom_at_max_roll,lf_diff_at_lift_onset,fbackminroll,fbackmaxroll,fbackmeanroll, fbackstdroll, roll_at_valley, roll_to_vall_time,timeNet_lf_initiate, timeNet_gf_initiate, gfLfCoupling])  
          
        if trialnum == starttrial:
            summary = trialsummary
        else:   
            summary = np.vstack([summary,trialsummary])

    #Save summary file
    #Create the subject's results directory, if it hasn't already
    if os.path.exists(os.path.join(DIR,'results',subject)) is False:
        os.makedirs(os.path.join(DIR,'results',subject))
        print('Created New Folder')
    #Save summary file    
    # # If created column(s), add name(s) here
    column_names = "totaltrialnum, com_cond, firsttriggertime, starttone, buttonuptime, onsettime, buttondowntime, Tx_onsetindx, onsety, Ty_onsetindx, Tz_onsetindx, rolloffset_onsetindx, absmaxroll, onsetofmaxroll, reachonset, reachphase, executeonset, executephase, Contact_COP_difference, gf_rate_thumb_max, gf_rate_pointer_max, lf_rate_thumb_max, lf_rate_pointer_max, COP_diff_at_lift_onset, COP_pointer_at_lift_onset, COP_thumb_at_lift_onset, gf_thumb_at_lift_onset, gf_pointer_at_lift_onset, lf_thumb_at_lift_onset, lf_pointer_at_lift_onset, load_phase, Mcom_at_lift_onset, Mcom_at_max_roll, lf_diff_at_lift_onset, fbackminroll, fbackmaxroll, fbackmeanroll, fbackstdroll, roll_at_valley, roll_to_vall_time, timeNet_lf_initiate, timeNet_gf_initiate, gfLfCoupling"
   
    np.savetxt(os.path.join(DIR,'results',subject,str(subject)+'summary_'+'starttrial'+str(trial_begin)+"_"+'.txt'), summary, fmt="%.3f", header=column_names, comments='')
    #np.savetxt(os.path.join(DIR,'results',subject,str(subject)+'summary_'+'starttrial'+str(trial_begin)+"_"+'.txt'),summary,fmt="%.3f",header=column_names, comments='')
