# adapted-aim2.py

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""11111111
Created on Thursday 17 July

@author: motorskill
CAS edits: aim 2 of R01
edits by jackson for VR protocol
"""

# imporing required python packages
import time
import serial
import viz, vizact, vizmat, viztask, viztracker, vizinfo, vizshape, vizmatplot
import numpy as np
import math, random
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

import warnings
warnings.filterwarnings("ignore")

import glob
import os.path
import nidaqmx
from nidaqmx import stream_readers
from nidaqmx.constants import AcquisitionType
from nidaqmx.constants import (TerminalConfiguration)

import fmri_trigger_test_2 as fmri_trigger
from fmri_trigger_test_2 import BUTTON_UP_EVENT

motion_freq = 150 # sampling frequency of motion data
force_rate = 500 # sampling frequency of force data
plot_redraw_rate = 15

import socket
import threading

class UnityBridge:
    def __init__(self, host='localhost', port=5555):
        self.sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        self.sock.bind((host, port))
        self.sock.listen(1)
        self.conn = None
        self.conn_lock = threading.Lock()

        # Start listener in background thread
        threading.Thread(target=self._accept_loop, daemon=True).start()
        print(f"[UnityBridge] Server started on {host}:{port}, waiting for Unity...")

    def _accept_loop(self):
        while True:
            conn, addr = self.sock.accept()
            print(f"[UnityBridge] Unity connected from {addr}")
            with self.conn_lock:
                self.conn = conn

    def wait_for(self, expected_msg):
        while True:
            with self.conn_lock:
                if self.conn:
                    try:
                        data = self.conn.recv(1024).decode('utf-8').strip()
                        print(f"[UnityBridge] Received: {data}")
                        if data == expected_msg:
                            return
                    except:
                        self.conn = None  # Reset connection on error
                time.sleep(0.1)

    def send(self, msg):
        with self.conn_lock:
            if self.conn:
                try:
                    self.conn.sendall((msg + "\n").encode('utf-8'))
                    print(f"[UnityBridge] Sent: {msg}")
                except:
                    print("[UnityBridge] Send failed. Resetting connection.")
                    self.conn = None

    def close(self):
        self.sock.close()
        with self.conn_lock:
            if self.conn:
                self.conn.close()

class MonkeyBars(object):
	X_OFFSET = 0
	Z_OFFSET = 0 
	Y_OFFSET = 0
	oldpos = [0,0,0]
	oldcount = 0
	force_bias_offset = [0,0,0,0,0,0,0,0,0,0,0,0]

	def __init__(self):
		"""Class that creates and enables interaction with
		objects needed for the flexion/extension + rotation
		study.
		"""
#		self.init_objects()
		
		# Calibration matrices for the force sensors.

		self.WM_32339 = np.array([[-0.263523377,	-0.053234635, -0.237917631,	6.801518886, 0.923811307, -7.589349989],\
                     [0.723376835, -8.101282297, -0.400267188, 3.880264378,	-0.344376527, 4.441885112],\
                         [15.94555466, -0.204333509, 15.9991545, 0.136270353, 16.88561928, 0.097599838],\
                             [0.004536927, -0.070985514, 0.170508814, 0.034014253, -0.179012639, 0.038804399],\
                                 [-0.193044668, 0.003981503, 0.097831088, -0.058571897,	0.098652848, 0.065761501],\
                                     [0.006408332, -0.07432435,	0.005129167, -0.068661303, 0.007118337,	-0.07885596]])

		self.WM_32340 = np.array([[-0.338236543, -0.022619792, 0.234400472,	6.974291267, 0.455601743, -7.333379685],\
                      [0.391093525, -8.026268253, -0.387468465,	4.009760144, -0.17152322, 4.242315924],\
                          [16.11414082, 0.039695741, 15.9906985, 0.078340837, 16.92229018, -0.015412285],\
                              [0.005045617,	-0.070050975, 0.171245309,	0.035192633, -0.176176261,	0.038228395],\
                                  [-0.1930104, -0.000170934, 0.10105247, -0.060532648,	0.09341937,	0.063977437],\
                                      [-0.000582218, -0.074812095,	0.003923866, -0.071382896,	0.004729957, -0.076053777]])

		self.init_objects()
		self.tracking_data = None
		self.tracking_data_two = None
		self.init_camera()
		self.time_capture = ()
		self.ati_voltage_data = []
		self.array_time_2 = list(range(4000))
		self.roll  = 0
		self.pitch = 0
		self.yaw = 0
		self.checkroll = 0
		self.checkrollvalues = False
		self.all_force_data = []
		self.plotready = False
		self.output = np.zeros([13, 4000])
		self.fdataready = False
		self.trial_number = None         
		
	def init_camera(self):
		"""Positions the camera such that all the experiment-
		relevant objects are visible
		"""
		viz.MainView.setPosition([0,0,-2.5])
		viz.MainView.setEuler([0,0,0],mode=viz.ABS_GLOBAL)   
		
	def init_objects(self):
		
		#ball (visible or not visible)
		self.ball = viz.addChild('white_ball.wrl',color=viz.PURPLE)
		self.ball.setScale([.5,.5,.5])
		self.ball.visible(viz.ON)
	
	def move(self,new_position):
		#move ball based on PPT position, centered based on offset calibration, and with gain of 5
		#New for simple visual: ball position should be based on constant moving in a circle; make ball visual
		#Disconnected to WorldViz
		self.ball.setPosition(new_position[0]*5,new_position[1]*5,new_position[2]*5)
		self.ballpos = self.ball.getPosition()
		if self.ballpos[0] == self.oldpos[0] and self.ballpos[1] == self.oldpos[1] and self.ballpos[2] == self.oldpos[2]:
			self.oldcount += 1
			if self.oldcount < 5:
				 'WARNING. NOT TRACKING'
		else:
			if self.oldcount > 0:
				 'TRACKING AGAIN'
			self.oldcount = 0
		self.oldpos = self.ballpos
		
	def tracker_rotate(self):
		"""This function collects the motion tracking data and is repeatedly called
        by the vizactontimer at the collection frequency."""
		#.01 on PPT = 1 cm
		
		#For your rigid-body tracker
		self.the_Eulers = tracker.getEuler() # Returns [yaw,pitch,roll] of node.
#		print('Live Roll (degrees):', self.the_Eulers[2])
		self.the_Positions = tracker.getPosition() # Returns [x,y,z] of node 
		#print('the_Positions',self.the_Positions)
		self.the_Positions[0] = self.the_Positions[0]+self.X_OFFSET
		self.the_Positions[1] = self.the_Positions[1]+self.Y_OFFSET
		self.the_Positions[2] = self.the_Positions[2]+self.Z_OFFSET
		
		self.move(self.the_Positions)
		
		# Also add the check roll here
		if self.checkrollvalues:		
			if experiment.initialroll < 0 and (self.the_Eulers[2])+ np.abs(experiment.initialroll) > 5:
				#flag bad roll 
				self.checkroll = 1
				#change background color on screen to red
				viz.clearcolor(1,0,0)
			
			if experiment.initialroll > 0 and np.abs(self.the_Eulers[2])- np.abs(experiment.initialroll) > 5:
				#flag bad roll 
				self.checkroll = 1
				#change background color on screen to red
				viz.clearcolor(1,0,0)
		
			else:
				#background color on screen black
				viz.clearcolor(0,0,0)
				
		self.the_Positions2 = tracker2.getPosition()
		self.the_Positions2[0] = self.the_Positions2[0]+self.X_OFFSET-self.the_Positions[0]
		self.the_Positions2[1] = self.the_Positions2[1]+self.Y_OFFSET-self.the_Positions[1]
		self.the_Positions2[2] = self.the_Positions2[2]+self.Z_OFFSET-self.the_Positions[2]
		
		self.the_Positions3 = tracker3.getPosition()
		self.the_Positions3[0] = self.the_Positions3[0]+self.X_OFFSET-self.the_Positions[0]
		self.the_Positions3[1] = self.the_Positions3[1]+self.Y_OFFSET-self.the_Positions[1]
		self.the_Positions3[2] = self.the_Positions3[2]+self.Z_OFFSET-self.the_Positions[2]
		
	def log_motion_data(self): 
		"""This function collects the motion tracking data and is repeatedly called
		by the vizactontimer at the collection frequency."""
		trackerdata = (self.the_Positions[0],self.the_Positions[1],self.the_Positions[2])
		trackerrotate = (self.the_Eulers[0],self.the_Eulers[1],self.the_Eulers[2])
		trackerdata2 = (self.the_Positions2[0],self.the_Positions2[1],self.the_Positions2[2])
		trackerdata3 = (self.the_Positions3[0],self.the_Positions3[1],self.the_Positions3[2])
		new_time_two = viz.tick()
		the_tracking_two = '\t'.join([str("%.3f"%new_time_two),
							str("%.3f %.3f %.3f"%trackerdata),str("%.3f %.3f %.3f"%trackerrotate),str("%.3f %.3f %.3f"%trackerdata2),str("%.3f %.3f %.3f"%trackerdata3)])
		write_data_two = the_tracking_two+'\n'
        # writing data to file opened in target order
		self.tracking_data_two.write(write_data_two)
		self.tracking_data_two.flush()
		
	def collect_force_data(self,trialnum):
		"""'This command (viz.director) executes the specified function in a separate thread, 
	thereby allowing you to perform operations that would otherwise halt rendering of the world or 
	the execution of other code within the script.' - Vizard'"""
		reader.read_many_sample(data = self.output, number_of_samples_per_channel = 4000)
		task.stop()
		self.all_force_data = []
		data_2 = self.output.T
		for row in data_2:
			result = row - self.force_bias_offset
			transducer1 = np.dot(self.WM_32339,(result[0:6].T))
			transducer2 = np.dot(self.WM_32340,(result[6:12].T))
			ar_volt = result[12:13].T
			transducer_combined = np.concatenate([transducer1,transducer2,ar_volt])
			transducer_combined_2 = transducer_combined.tolist()
			self.all_force_data.append(transducer_combined_2)
		self.force_array = np.array(self.all_force_data)
		self.force_array_transpose = self.force_array.T
		self.output = np.zeros([13, 4000])
		self.fdataready = True
		
	def task_plotting(self):
		"""Task that updates plot with latest data and redraws texture. Called in seperate schedule in main function."""
		while True:
			yield viztask.waitTime(1.0 / plot_redraw_rate)
			if self.plotready:
				array_fy1 = self.force_array_transpose[1]
				array_fz1 = self.force_array_transpose[2]
				array_fy2 = self.force_array_transpose[7]
				array_fz2 = self.force_array_transpose[8]
				line_fy1.set_data(self.array_time_2,array_fy1)
				line_fz1.set_data(self.array_time_2,array_fz1)
				line_fy2.set_data(self.array_time_2,array_fy2)
				line_fz2.set_data(self.array_time_2,array_fz2)
				ax_forces_grip.set_title('Trial ' + self.trial_number)                 
				# redraws the figure 
				yield fig_texture.waitRedraw()
				self.plotready = False


        
class find_Tracker(viztask.Condition ):
	#This class will return when the participant has moved within the range of the correct position. 
	#Origin is the point of interest
	#Threshold is the allowance for how close the subject can be
	#Need to play around with this so it's accurate but doesn't miss positions
	#CheckMove = find_Tracker(bars,.02,[0,0.03])
	#movetoStart = yield viztask.waitAll([CheckMove])
	def __init__(self,bars,threshold,origin):
		self.bars = bars
		self.__thres = threshold
		self.__yori = origin

	def update(self):
		pos = bars.the_Positions
		print("y position =", pos[1])
		if pos[1] < self.__thres:
			print("threshold in y direction has been met")
			return True
		return False
		
class Experiment():	

	def __init__(self,bars):

		self.bars= bars
		self.tracking_data = None
		self.tracking_data_two = None
		self.logfile2 = None
		self.firstrealtrigger = None
		self.task = None  
		self.id = None
		self.trial_begin = None            
		self.randorder = []
		self.trialnum = [] 
		self.itivalue = []
		self.COM_cond = []  
		self.repcondition=[]
		#self.LED_cond = []
		#self.pairType = []
		self.experiment_data = None
		self.target_threshold = 0.01 # 0.01 is 1cm 
		self.unity_bridge = None

        # Trialwise data
		self.tsyncsig = -1
		self.starttone = -1
		self.endtone = -1
		self.tdown = -1
		self.treallydown = -1
		self.eventnum = -1
		self.triggered = -1
		self.button_up = -1
		self.button_down = -1
		self.startiti = -1
		self.enditi = -1

	def participantInfo(self):
		#add the vizinfo box
		infoBox = vizinfo.add('')
		infoBox.scale(1,2)
		infoBox.translate(0.4,0.8)
		infoBox.title('Participant Info - Main Experiment')

		participantID = infoBox.add(viz.TEXTBOX,'Participant ID')
		start_trial = infoBox.add(viz.TEXTBOX, 'Starting Trial')
		condition_type = infoBox.add(viz.TEXTBOX, 'Condition')
		runButton = infoBox.add(viz.BUTTON,'Run')
		yield viztask.waitButtonUp(runButton)

		self.id = participantID.get()
		self.trial_begin = start_trial.get()
		self.condition_type = condition_type.get()
		infoBox.remove()
	
	def get_conditions(self):        
		filename = glob.glob(path+'Trials/repVR/'+str(self.id)+"_"+str(self.condition_type)+'*.txt') 
		print(filename)		
		file = open(filename[0],'r')
		lines = file.readlines()[1:]

		for line in lines:                 
			s=line.split(",")
			print(s)
			self.trialnum.append(float(s[3]))
			self.COM_cond.append(float(s[0]))
			self.itivalue.append(float(s[2]))
			self.repcondition.append(float(s[1]))			
			
			
		trials = len(self.trialnum)
		startrow = int(self.trial_begin) - 1 
	
		self.randorder = list(range(startrow,trials))
		file.close() #edit the next bit definin the path
		if not os.path.exists('Data\\repVR\\S'+str(self.id)+"_"+str(self.condition_type)):
			os.makedirs('Data\\repVR\\S'+str(self.id)+"_"+str(self.condition_type))
		exp_file = 'Data\\repVR\\S'+str(self.id)+'_'+str(self.condition_type)+'\\S'+str(self.id)+"_"+str(self.condition_type)+'_experiment_data'+str(RANDOM)+'_'+str(self.trial_begin)+'.txt'
		self.experiment_data = open(exp_file,'a')
		
		#if self.LED_cond[0] == 1 or self.LED_cond[0] == 2:
		#	ser.write(b'H')           
		#elif self.LED_cond[0] == 3 or self.LED_cond[0] == 4:
		#	ser.write(b'L')
			
		
		#if self.LED_cond[0] == 1 or self.LED_cond[0] == 3:
		#	ser.write(b'O')
		#elif self.LED_cond[0] == 2 or self.LED_cond[0] == 4:
		#	ser.write(b'F')

        # Checking the center of mass condition for the first trial. 
        # Send a signal to the Arduino to send the wight to the or to the left with the pneumatic system
        			

		print('COM_cond: ',self.COM_cond[0])
		
	def set_origin(self):
		print("Waiting for mousedown to calibrate origin")
		yield viztask.waitMouseDown(None)

		self.bars.X_OFFSET = -bars.the_Positions[0]
		self.bars.Y_OFFSET = -bars.the_Positions[1]
		self.bars.Z_OFFSET = -bars.the_Positions[2]

		print("waiting for another mouse down")
		yield viztask.waitMouseDown(None)
		print("continuing")
	
		#save offset values in .txt file #change the savename here!!!
		alloffsets = np.array([bars.X_OFFSET,bars.Y_OFFSET,bars.Z_OFFSET])
		np.savetxt('Data\\repVR\\S'+str(self.id)+"_"+str(self.condition_type)+\
		           '\\S'+str(self.id)+'_'+ \
				    '_'+'se_'+ \
					'_'+ \
					'offset_data' +str(RANDOM)+'_'+str(self.trial_begin)+ 
					'.txt',alloffsets,fmt="%.3f\t")

	def zero_out_force(self):
		"""Collects the bias offset voltage values. Averages 0.5 second of data for each channel."""
		print("Waiting for mousedown to Zero Out Force Sensors")
		yield viztask.waitMouseDown(None)
		bias_off_force_chunk = np.zeros([13, 500])
		reader.read_many_sample(data = bias_off_force_chunk, number_of_samples_per_channel = 500)
		task.stop()
		self.bars.force_bias_offset= bias_off_force_chunk.mean(axis = 1)
		print('force offset:',self.bars.force_bias_offset)
		print("continuing")


	def intertrial(self,i,iti):
		self.startiti = viz.tick()
		#CHANGE
		#wait for iti - 1 second
		yield viztask.waitTime(iti-1) # wait intertrial time - 1 
		#self.light_intertrial(i) # reset the light - turn on or off
        
		# Setting the Arduino to DAQ pin to low
		ser.write(b'N')
        
        # Cehcking the COM condition for the next trial

			
	#	print('COM_cond: ',self.COM_cond[i+1])
		
		#stop logging data
		timer_log_motion.setEnabled(viz.OFF)
		
		# MIGHT NEED TO ADD MORE TIME HERE DEPEDNING ON THE TIME IT TAKES TO TURN THE TABLE
		yield viztask.waitTime(1) # wait for 1 second
		self.enditi = viz.tick()

	def trial(self,i,trialnum):
		#1. tone indicating start
		#2. Time to do task
		#3. tone to put down
		#play tone
		bars.fdataready = False
		if self.COM_cond[i] ==1:  
			ser.write(b'P') # This sends a signal that sends the weight to the left
			#P AND S SWAPPED FOR TROUBLESHOOTING
			viz.playSound("resources/left_side.wav") # Audio that indicates the weight will be on the left side for the next trial        
		elif self.COM_cond[i] == 2:  
			ser.write(b'S') # This sends a signal that sends the weight to the 
			#P AND S SWAPPED FOR TROUBLESHOOTING
			viz.playSound("resources/right_side.wav")
		yield viztask.waitTime(1)
		viz.director(self.bars.collect_force_data,trialnum)   # Opens seperate thread for transducer voltage data log
		timer_log_motion.setEnabled(viz.ON) # turns on timer for motion data log
		
		yield viztask.waitTime(0.25) # ADDED this wait time so that event occurs after force starts logging 
		ser.write(b'T')
		self.tsyncsig = viz.tick()
#		
		
		bars.checkrollvalues = False
		
		#synch start to TR
		print("Waiting for a trigger event")
		yield viztask.waitEvent(fmri_trigger.TRIGGER_EVENT)
		#tone to start
		viz.playSound("resources/beep-louder.wav") # pick up tone
		self.initialroll = bars.the_Eulers[2]
		print(self.initialroll)
		print("Got a trigger")
		self.eventnum = trigger.event_count  #TR num of start of trial
		self.triggered = trigger.triggertime  #timestamp of start of trial
		self.starttone = viz.tick()

		print("waiting for button up L271")
		yield viztask.waitEvent(BUTTON_UP_EVENT) # waits for button box button to be released
		# ser.write(b'L') # Turn off the light at button-up
		print("button got released")
		self.button_up = viz.tick()
		bars.checkrollvalues = True
		
		# 1 second to grasp object
		#yield viztask.waitTime(1)
		#viz.playSound("resources/beep-01.wav")

		#2.5 seconds to do rest of task
		print("waiting 2.5 sec")
		yield viztask.waitTime(2.5) # Hold on object for 2.5 seconds
		print("2.5 sec over.")
		print('roll value:',bars.checkroll)
		yield bars.checkroll

		#tone to put down
		viz.playSound("resources/beep-07.wav") # put down tone
		self.endtone = viz.tick()
		bars.checkrollvalues = False
		
		#Check to see if T y position really is down.
		CheckMove = find_Tracker(bars,self.target_threshold,0)  #0.01 is your threshold, 0 is your target; 
		print("Waiting for tracker to get less than %.03f meters from the table" % self.target_threshold)
		print("OR a mousedown")
		yield viztask.waitAny([CheckMove, viztask.waitMouseDown(None)]) #wait for object to get close to table or mousedown
		self.treallydown = viz.tick()
		print("They put the object down")
		
		print("waiting for fmri_trigger.BUTTON_DOWN_EVENT", self.endtone)
		yield viztask.waitEvent(fmri_trigger.BUTTON_DOWN_EVENT) # waits for button down event
		self.button_down = viz.tick() 
		print("got it", self.button_down)
		
		if bars.checkroll==1: # play a specific warning sound if the roll during the trial was too much
			viz.playSound("resources/beep9.wav")
			bars.checkroll = 0 
			print('roll value:',bars.checkroll)

	def write_data(self,i,trialnum):
         while True:
              if bars.fdataready == False:
                   yield viztask.waitTime(.2)
              elif bars.fdataready == True:   
                   bars.trial_number = str(trialnum)                    
                   bars.plotready = True
                   OnsetInfo = '\t'.join([str(i+1),str(self.eventnum)])
                   TrialInfo = '\t'.join([str("%.3f"%self.itivalue[i]),str(self.COM_cond[i]),str(self.repcondition[i])])
				   #TrialInfo = '\t'.join([str("%.3f"%self.itivalue[i]),str(self.COM_cond[i]),str(self.LED_cond[i]),str(self.pairType[i])])
                   Timestamps = '\t'.join([str("%.3f"%self.tsyncsig),str("%.3f"%self.firstrealtrigger),str("%.3f"%self.triggered),str("%.3f"%self.starttone),str("%.3f"%self.endtone),str("%.3f"%self.treallydown),str("%.3f"%self.startiti),str("%.3f"%self.enditi),str("%.4f"%self.button_up),str("%.3f"%self.button_down)])
    	            #took out str("%.3f"%confirmendpoint)
    		
                   # writes out data onto experiment data text file (opned in get_conditions function)
                   # trial, number, Trigger event number at start, itivalue, light_delay, first trigger in experiment, start of trial time, start of force collection time, endtone time, really_down time, start of intertrial time, end of intertrial time, end of intertrial time, button_up time, button_down time
                   self.experiment_data.write(OnsetInfo+'\t'+TrialInfo+'\t'+Timestamps+'\n')
                   self.experiment_data.flush()
            
					# saves a text file for all the appended trigger pulse events (event # and timestamp) at this point. Appended in trigger module.
                   np.savetxt('Data\\repVR\\S'+str(self.id)+"_"+str(self.condition_type)+'\\S'+str(self.id)+"_"+str(self.condition_type)+'_'+str(int(trialnum))+'_'+'trigger_data' +str(RANDOM)+'_'+str(self.trial_begin)+'.txt',np.vstack((trigger.events,trigger.rt)).T,fmt="%.3f\t") #add condition after this
                   # saves a text file for all the appended trigger pulse events (event # and timestamp) at this point. Appended in trigger module.
                   np.savetxt('Data\\repVR\\S'+str(self.id)+"_"+str(self.condition_type)+'\\'+'S'+str(self.id)+"_"+str(self.condition_type)+'_'+str(self.trialnum[i])+'_' +'force_data'+str(RANDOM)+'_'+str(self.trial_begin)+'.txt', self.bars.force_array, delimiter=",")
                   self.bars.force_array = []
                   break
		
	def target_order(self):
	
		print('path=',path)		
		yield self.participantInfo() # Opens graphical user interface
		yield self.get_conditions()  # Gets conditions from text file (intertrial time - delay for light)
		yield self.zero_out_force()  # Gets bias offset for voltage
		yield self.set_origin()      # sets the bias offset for markers 
		yield viztask.waitTime(2)    # Wait two seconds 

		self.unity_bridge.send("begin_experiment")
		
		trigger.start()              # Starting trigger module
		print('TURN ON SCANNER')
		TRstowait = 15 #15*.4 = 6 s
		print('wait 6 s (15 TRs) at beginning -- should not begin counting unless scanner is on!')
		trialstart = 0
		while trialstart < TRstowait: # Buffer of 15 triggers 
			yield viztask.waitEvent(fmri_trigger.TRIGGER_EVENT)
			if trialstart == 0:
				self.firstrealtrigger = viz.tick()  #time in program when the first scanner trigger came in
			trialstart += 1
			print(trialstart)
		
		yield unity_bridge.wait_for("unity_ready")  # wait at the start
		
		for i in self.randorder: # for each trial 
			if i in [40]: # After 10 switches pause for thirty seconds, after each set of 20 trials pause for 1 minute.  
				print("Break")
				yield viztask.waitMouseDown(None)
            
            # Opens up file to store motion tracking data
			logfile2 = path+ 'Data\\repVR\\S'+str(self.id)+"_"+str(self.condition_type)+'\\'+'S'+str(self.id)+"_"+str(self.condition_type)+'_'+str(self.trialnum[i])+'_motion_data'+str(RANDOM)+'_'+str(self.trial_begin)+'.txt'
			self.bars.tracking_data_two = open(logfile2,'a')
			yield self.trial(i,self.trialnum[i]) # Goes through trial function
			yield self.intertrial(i,self.itivalue[i]) #Goes through intertrial
			yield self.write_data(i,self.trialnum[i]) 
		endDisplay = viz.addText( 'Run Completed', viz.SCREEN)
		endDisplay.alignment(viz.ALIGN_CENTER)
		endDisplay.visible(1)
		endDisplay.visible(1)
		endDisplay.setPosition(.5,.65)
		endDisplay.fontSize(45)
		print('wait 10 s (25 TRs) at end')
		TRstowait = 25
		trialstart = 0
		while trialstart < TRstowait: #buffer of 25 trigger pulses at the end
			yield viztask.waitEvent(fmri_trigger.TRIGGER_EVENT)
			trialstart += 1
			print(trialstart)
		unity_bridge.send("quit")  # at the end
		unity_bridge.close()
		endDisplay.visible(0)
		viz.quit()

if __name__=="__main__":
    """ This is the main function of the script. """
    unity_bridge = UnityBridge()
	
    viz.go() # starts Vizard rendering environment
    viz.fov(25) # Setting for window that it opens
    viz.vsync(0) # turns off vertical synchronizations (Does not need to be matched with monitor)
    viz.setOption('viz.max_frame_rate',150) # global variable for maximum frame rate
    viz.window.setPosition(620,80)
    viz.window.setSize(1300,660)
    viz.window.setFloat()
    viz.window.setName( 'Inteface and Feedback Window')

    window = viz.MainWindow
    MainView = viz.addView()
    window.setView(MainView)
    MainView.setPosition([1.2,0,-6])

    SecondWindow = viz.addWindow(pos=(0.5,1),size=(0.5,1))
    SecondWindow.visible(0,viz.SCREEN)

    SecondView = viz.addView()
    SecondWindow.setView(SecondView)
    SecondView.setPosition([0,5,-1])
    SecondWindow.fov(60)

    vizinfo.InfoPanel('Object Motion Feedback', window = viz.MainWindow)
    vizinfo.InfoPanel('Digit Forces Feedback', window = SecondWindow)

    path = 'C:\\Users\\WorldViz\\Documents\\Python Scripts\\' # path for files
    RANDOM = random.randint(1001,9999) # Rndom number used in text file names for each run
    
    viz.playSound("resources/beep-07.wav",viz.SOUND_PRELOAD) #Preloading wav files for later
    viz.playSound("resources/beep-01.wav",viz.SOUND_PRELOAD)
    viz.playSound("resources/beep-louder.wav",viz.SOUND_PRELOAD)
    viz.playSound("resources/left_side.wav",viz.SOUND_PRELOAD)
    viz.playSound("resources/right_side.wav",viz.SOUND_PRELOAD)
    viz.playSound("resources/beep9.wav",viz.SOUND_PRELOAD)

    # opening a serial port connection to microcontroller
    # Arguments = communication port,rate, timelimit for writing
    ser = serial.Serial('COM3',115200,writeTimeout=0) 
    time.sleep(2) # allow time for port connection to be established

	#call Monkey Bars Class
    bars = MonkeyBars()
    experiment = Experiment(bars=bars)
    experiment.initialroll = 0
    experiment.unity_bridge = unity_bridge

	#Turn off mouse navigation for vizard window
    viz.mouse.setOverride()
	
    # "A task is a collection of one or more virtual channels with timing, triggering, and other properties." - NI
    task = nidaqmx.Task()
    task.ai_channels.add_ai_voltage_chan("Dev1/ai0", terminal_config=TerminalConfiguration.DIFFERENTIAL) #Defining analog inputs at input channels
    task.ai_channels.add_ai_voltage_chan("Dev1/ai1", terminal_config=TerminalConfiguration.DIFFERENTIAL)
    task.ai_channels.add_ai_voltage_chan("Dev1/ai2", terminal_config=TerminalConfiguration.DIFFERENTIAL)
    task.ai_channels.add_ai_voltage_chan("Dev1/ai3", terminal_config=TerminalConfiguration.DIFFERENTIAL)
    task.ai_channels.add_ai_voltage_chan("Dev1/ai4", terminal_config=TerminalConfiguration.DIFFERENTIAL)
    task.ai_channels.add_ai_voltage_chan("Dev1/ai5", terminal_config=TerminalConfiguration.DIFFERENTIAL)
    task.ai_channels.add_ai_voltage_chan("Dev1/ai16", terminal_config=TerminalConfiguration.DIFFERENTIAL)
    task.ai_channels.add_ai_voltage_chan("Dev1/ai17", terminal_config=TerminalConfiguration.DIFFERENTIAL)
    task.ai_channels.add_ai_voltage_chan("Dev1/ai18", terminal_config=TerminalConfiguration.DIFFERENTIAL)
    task.ai_channels.add_ai_voltage_chan("Dev1/ai19", terminal_config=TerminalConfiguration.DIFFERENTIAL)
    task.ai_channels.add_ai_voltage_chan("Dev1/ai20", terminal_config=TerminalConfiguration.DIFFERENTIAL)
    task.ai_channels.add_ai_voltage_chan("Dev1/ai21", terminal_config=TerminalConfiguration.DIFFERENTIAL)
    task.ai_channels.add_ai_voltage_chan("Dev1/ai23", terminal_config=TerminalConfiguration.RSE)	
    task.timing.cfg_samp_clk_timing(rate = force_rate,sample_mode = AcquisitionType.CONTINUOUS) # Sampling frequency set
    reader = stream_readers.AnalogMultiChannelReader(task.in_stream)
    task.start()

    array_time, array_fy1,array_fz1,array_fy2,array_fz2  = [],[],[],[],[]
    fig = plt.figure()

    ax_forces_grip = fig.add_subplot(211)#create axes
    ax_forces_grip.axis([0,4000,-25,5])  #set axis limits UPDATE ME AFTER TROUBLESHOOT
    ax_forces_grip.set_xticklabels([])	
	

    ax_forces_lift = fig.add_subplot(212)#create axes
    ax_forces_lift.axis([0,4000,0,12])  #set axis limits
    ax_forces_lift.set_xticklabels([])	

    line_fy1, = ax_forces_lift.plot(array_time,array_fy1,'tab:blue', label = 'Sensor1')
    line_fz1, = ax_forces_grip.plot(array_time,array_fz1,'tab:blue', label = 'Sensor1')
    line_fy2, = ax_forces_lift.plot(array_time,array_fy2,'tab:orange', label = 'Sensor2')
    line_fz2, = ax_forces_grip.plot(array_time,array_fz2,'tab:orange', label = 'Sensor2')

    ax_forces_grip.set_ylabel('Grip (N)')
    ax_forces_lift.set_ylabel('Lift (N)')
    ax_forces_grip.legend(frameon = False);
    ax_forces_lift.legend(frameon = False);

    fig_texture = vizmatplot.Texture(fig)

	# Create quad to render plot texture
    quad = viz.addTexQuad(texture=fig_texture, pos=(0,5,0))
    quad.alpha(0.8)

	# Setting up motion tracking capabilities 
    vrpn = viz.add('vrpn7.dle') # Connecting to PPT Data within Vizard 
    tracker = vrpn.addTracker('PPT0@localhost',0) # rigid_body tracker (composed of markers 0 through 2)
    tracker2 = vrpn.addTracker('PPT0@localhost',3)  # position tracker for marker 3
    tracker3 = vrpn.addTracker('PPT0@localhost',4)  # position tracker for marker 4

	# Setting up python module to be used for collecting trigger pulses ...
	# ... and creating events for button box pushes.

    fmrion = False #CHANGE TO True WHEN USING SCANNER!!!!!!!
    using_celeritas = True

    if fmrion:
        trigger = fmri_trigger.FmriTrigger()
    else:
        if using_celeritas:
            trigger = fmri_trigger.FakeFmriTrigger()
        else:
            trigger = fmri_trigger.TotallyFakeFmriTriggerWithMouse()

	# Motion tracking timers for data collection and data logging
	# Arguments = rate (seconds), # of times to repeat, function being called
    timer_motion = vizact.ontimer2((1/motion_freq), viz.FOREVER,bars.tracker_rotate)
    timer_motion.setEnabled(viz.ON)
    timer_log_motion = vizact.ontimer2((1/motion_freq),viz.FOREVER,bars.log_motion_data)
    timer_log_motion.setEnabled(viz.OFF)

    # Schedules task and removes it when finished.
    # Controls data collection, sounds, and lights during run of trials
    viztask.schedule(bars.task_plotting())
	
    # Schedules task and removes it when finished.
    # Controls force plots during experiment
    viztask.schedule(experiment.target_order())
