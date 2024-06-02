It's a intigrated navigation proj used by st-ekf.

estimator.cpp is for loose couple and tightly_estimator.cc is for tightly couple. 


# Build and implement
cd Build 
make
cmake
#Result
loose couple result is saved at zuhe_pos_vel_yaw1.txt, tightly couple result is saved at tight_estimator_result.txt.
#plot
python plot_result.py
python plot_tightly_result.py
#for your dataset
make sure the dataset path is right.
