import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap


DATA_URL = 'http://earthquake.usgs.gov/earthquakes/feed/v1.0/summary/4.5_month.csv'
print("Downloading", DATA_URL)
df = pd.read_csv(DATA_URL)

fig, ax = plt.subplots()
earth = Basemap(ax=ax)
earth.drawcoastlines(color='#556655', linewidth=0.5)
ax.scatter(df['longitude'], df['latitude'], df['mag'] ** 2, 
           c='red', alpha=0.5, zorder=10)
ax.set_xlabel("This month's 4.5M+ earthquakes")
fig.savefig('usgs-monthly-4.5M.png')
