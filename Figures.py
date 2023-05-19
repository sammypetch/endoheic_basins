
# Flux-inferred storage figure

plt.figure(figsize=(10,5))
plt.plot(x,som_stor2020[0:216],c='orange',label='GRACE')
plt.plot(x,FIS_PEobsx, c='purple',label='FIS')
plt.legend(fontsize=12)
plt.ylabel('Total water storage anomaly (cm)', fontsize=13)
plt.xlim(2002,2020)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.ylim(-55,12)
plt.xlabel('Time', fontsize=13)
plt.tight_layout()
plt.show()

# mean seasonal storage cycle

plt.figure()
plt.ylabel('cm', fontsize=12)
plt.xlabel('Month', fontsize=12)
plt.plot(x12,mscFISPopt2, c = 'b', label = 'FIS (P) ')
plt.plot(x12,mscFISPobs2, c = 'b', ls='dashed')
plt.plot(x12,mscFISLEopt2, c = 'red', label = 'FIS (- E)')
plt.plot(x12,mscFISLEobs2, c = 'red', ls='dashed')
plt.plot(x12,mscFISPEopt2 , c = 'purple', label = 'FIS (P - E) ')
plt.plot(x12,mscFISPEobs2, c = 'purple', ls='dashed')
plt.plot(x12,mscgrace, c= 'orange', label = 'GRACE')
plt.legend(fontsize=11.5)
plt.xticks(x12,months, fontsize=12)
plt.yticks(fontsize=12)
plt.ylim(-5,5)
plt.xlim(0,11)
plt.show()

# mean seasonal flux cycle

plt.figure()
plt.ylabel('mm day$^{-1}$ ', fontsize=12)
plt.xlabel('Month', fontsize=12)
plt.plot(x12,mmPopt, c = 'b', label = 'P')
plt.plot(x12,mmPobs, c = 'b', ls='dashed')
plt.plot(x12,mmLEopt, c = 'red', label = 'E')
plt.plot(x12,mmLEobs, c = 'red', ls='dashed')
plt.plot(x12,mmPEopt , c = 'purple', label = 'P - E')
plt.plot(x12,mmPEobs, c = 'purple', ls='dashed')
plt.plot(x12,mmdSGRACE, c= 'orange', label = 'GRACE dS')
plt.legend( fontsize =12)
plt.xticks(x12,months, fontsize=12)
plt.yticks(fontsize=12)
plt.xlim(0,11)
plt.show()



# inter-annual variability plot

plt.figure(figsize =(10,5))
plt.plot(x, FIS_PoptDS, c ='b', label = 'FIS (P)')
plt.plot(x, FIS_LEoptDS, c ='r', label = 'FIS (-E)')
plt.plot(x, FIS_PobsDS, c ='b', ls='dashed')
plt.plot(x, FIS_LEobsDS, c ='r', ls='dashed')
plt.plot(x, FISgDS, c ='orange', label = 'GRACE')
plt.xlim(2002,2020)
plt.ylabel('Total water storage anomaly (cm)', fontsize=12)
plt.xlabel('Time', fontsize=12)
plt.legend(fontsize=13)
plt.yticks(fontsize=12)
plt.xticks(fontsize=12)
plt.xlim(2002,2020)
plt.show()


# Energy storage figure

plt.figure(figsize=(10,5))
plt.plot(x, energy_Stor_obs, ls = 'dashed', label = 'Observations')
plt.plot(x, energy_Stor_opt, c = 'green', label = 'After yearly budget closure')
plt.legend(fontsize=12)
plt.xlim(2002,2020)
plt.ylim(-10e8,7e8) # setaccording to basin
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.ylabel('Total energy storage anomaly (Jm$^{-2}$)', fontsize=12)
plt.xlabel('Time', fontsize=12)
plt.show()



# Trend Figure

fig , ax1 = plt.subplots(figsize=(10,5))
plt.xlim(2002,2019)
ax1.set_xlabel('Time', fontsize=12)
plt.yticks(fontsize=12)
plt.xticks(fontsize=12)
ax2 = ax1.twinx()
ax1.plot(Xy, Pcmyear, c = 'b')
ax1.plot(Xy, Ptrend, c = 'b', ls='dashed', label = 'Total storage change from P trend = ')
ax1.scatter(Xy, Pcmyear, c = 'b')
ax1.set_ylabel('Annual mean P (cm year$^{-1}$)', fontsize =12, c='b' )
ax1.set_ylim(17,25)
ax1.tick_params('y',colors='b')
plt.yticks( fontsize = 12)
ax2.plot(Xy, dScmyear, c = 'orange')
ax2.scatter(Xy, dScmyear, c= 'orange')
ax2.plot(Xy, dStrend, c= 'orange', ls = 'dashed', label='Total storage change from dS trend = ')
ax2.set_ylabel('Annual mean dS (cm year$^{-1}$)',fontsize= 12, c='orange')
ax2.set_ylim(-4,4)
ax2.tick_params('y',colors='orange')
plt.xlabel('Time', fontsize =12)

handles1, labels1 = ax1.get_legend_handles_labels()
handles2, labels2 = ax2.get_legend_handles_labels()
# Combine the handles and labels
handles = handles1 + handles2
labels = labels1 + labels2
# Create the combined legend
plt.legend(handles, labels, loc='upper left', fontsize=13)
plt.tight_layout()
plt.savefig('as1_TRENDnew.jpg',dpi=300)
plt.show()
