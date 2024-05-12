import matplotlib.pyplot as plt
import numpy as np

with open("./tight_estimator_result.txt") as f:
    datas = []
    for line in f.readlines():
        data = []
        for pos_vel in line.split():
            data.append(float(pos_vel))
        datas.append(data)

with open("./gps_data.txt") as f:
    gnssdatas = []
    for line in f.readlines():
        gnssdata = []
        for pos_vel in line.split():
            gnssdata.append(float(pos_vel))
        gnssdatas.append(gnssdata)

error_long = np.subtract([pos_vel[1] for pos_vel in gnssdatas],[pos_vel[1] for pos_vel in datas])
error_lat = np.subtract([pos_vel[0] for pos_vel in gnssdatas],[pos_vel[0] for pos_vel in datas])

fig1, ax1 = plt.subplots()
ax1.plot([pose[7] for pose in datas], [pose[6] for pose in datas])
ax1.set_title('trajectory')
plt.xlabel('x')
plt.ylabel('y')
plt.show()

plt.xlim(-1e5, 1e5)
sub_fig1 = plt.subplot(2, 1, 1)
sub_fig1.plot([time for time in range(len(error_long))],[error for error in error_long])
sub_fig1.set_title('error_long')
plt.xlabel('time')
plt.ylabel('error deg')
sub_fig1.ticklabel_format(style='sci', axis='both', scilimits=(0, 0))

sub_fig2 = plt.subplot(2, 1, 2)
sub_fig2.plot([time for time in range(len(error_lat))],[error for error in error_lat])
sub_fig2.set_title('error_lat')
plt.xlabel('time')
plt.ylabel('error deg')
sub_fig2.ticklabel_format(style='sci', axis='both', scilimits=(0, 0))
fig1.subplots_adjust(hspace=5)
plt.show()