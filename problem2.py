import pygmt
import matplotlib.pyplot as plt
import rioxarray

dataarray =  rioxarray.open_rasterio('./dataset/displacement_dataset/ldem_4.tif')
tif_shape = dataarray.shape
dataarray = dataarray.isel(band=0)
import xarray as xr
import numpy as np 

height, width = dataarray.shape

xprecision = 360 / width 
yprecision = 180 / height

new_grid = xr.DataArray(
    data=dataarray,
    coords=dict(
        y=np.linspace(start=-90, stop=90, num=height),
        x=np.linspace(start=-180, stop=180, num=width),
    ),
    dims=("y", "x"),
)

color_dataarray =  rioxarray.open_rasterio('./dataset/color_dataset/lroc_color_poles_4k.tif')
color_tif_shape = color_dataarray.shape
color_dataarray = color_dataarray.isel(band=0)
import xarray as xr
import numpy as np 

color_height, color_width = color_dataarray.shape

color_grid = xr.DataArray(
    data=color_dataarray,
    coords=dict(
        y=np.linspace(start=-90, stop=90, num=color_height),
        x=np.linspace(start=-180, stop=180, num=color_width),
    ),
    dims=("y", "x"),
)

print(color_grid)


region = "90/110/5/35"
minlong, maxlong, minlat, maxlat = region.split('/')
minlong, maxlong, minlat, maxlat = float(minlong), float(maxlong), float(minlat), float(maxlat)
perspective = [150, 30]
selected_long = 0
selected_lat = 0
selected_point_height = new_grid[selected_long, selected_lat].data
contour_height = 1
projection = ''
cmap = 'gray'

whole_fig = pygmt.Figure()
whole_fig.grdimage(grid=color_grid, cmap="gray", projection="Q12C")

whole_fig.savefig('wholefig.png')

print('saved')

zoomed_fig = pygmt.Figure()
zoomed_fig.grdview(
    grid=new_grid,
    # Sets the view azimuth as 130 degrees, and the view elevation as 30
    # degrees
    perspective=perspective,
    # xmin/xmax/ymin/ymax
    region=region,
    # Sets the x- and y-axis labels, and annotates the west, south, and east
    # axes
    frame=["xa", "ya", "WSnE"],
    # Sets a Mercator projection on a 15-centimeter figure
    projection="Q12C",
    # Sets the height of the three-dimensional relief at 1.5 centimeters
    zsize="4.5c",
    surftype="s",
    cmap=cmap
)

zoomed_fig.colorbar(frame=["a1", "x+lElevation", "y+lm"])
zoomed_fig.savefig('tempfig.png')

zoomed_flat_fig = pygmt.Figure()
zoomed_flat_fig.grdimage(grid=new_grid, cmap=cmap, projection="Q12C", region=region)

zoomed_flat_fig.savefig('tempflatfig.png')
import cv2 

img = cv2.imread('wholefig.png')
resized = cv2.resize(img, (600,300))
cv2.imwrite('wholefig.png', resized)

img = cv2.imread('tempfig.png')
resized = cv2.resize(img, (600,300))
cv2.imwrite('tempfig.png', resized)

img = cv2.imread('tempflatfig.png')
resized = cv2.resize(img, (600,300))
cv2.imwrite('tempflatfig.png', resized)

def change_view():
    import os
    os.remove('tempfig.png')
    azimuth, elevation = perspective[0], perspective[1]
    azimuth = int(azimuth)
    elevation = int(elevation)
    zoomed_fig = pygmt.Figure()
    zoomed_fig.grdview(
        grid=new_grid,
        # Sets the view azimuth as 130 degrees, and the view elevation as 30
        # degrees
        perspective=[azimuth, elevation],
        region=region,
        # Sets the x- and y-axis labels, and annotates the west, south, and east
        # axes
        frame=["xa", "ya", "WSnE"],
        # Sets a Mercator projection on a 15-centimeter figure
        projection="Q12C",
        # Sets the height of the three-dimensional relief at 1.5 centimeters
        zsize="4.5c",
        surftype="s",
        cmap=cmap
    )
    
    zoomed_fig.colorbar(frame=["a1", "x+lElevation", "y+lm"])

    zoomed_fig.savefig('tempfig.png')

    img = cv2.imread('tempfig.png')
    resized = cv2.resize(img, (600,300))
    cv2.imwrite('tempfig.png', resized)

    print('done')

def change_zoomed_region(region):
    import os
    os.remove('tempfig.png')
    
    zoomed_fig = pygmt.Figure()
    zoomed_fig.grdview(
        grid=new_grid,
        # Sets the view azimuth as 130 degrees, and the view elevation as 30
        # degrees
        perspective=perspective,
        region=region,
        # Sets the x- and y-axis labels, and annotates the west, south, and east
        # axes
        frame=["xa", "ya", "WSnE"],
        # Sets a Mercator projection on a 15-centimeter figure
        projection="Q12C",
        # Sets the height of the three-dimensional relief at 1.5 centimeters
        zsize="4.5c",
        surftype="s",
        cmap=cmap
    )
    zoomed_fig.colorbar(frame=["a1", "x+lElevation", "y+lm"])
    zoomed_flat_fig = pygmt.Figure()
    zoomed_flat_fig.grdimage(grid=new_grid, cmap=cmap, projection="Q12C", region=region)

    zoomed_flat_fig.savefig('tempflatfig.png')
    img = cv2.imread('tempflatfig.png')
    resized = cv2.resize(img, (600,300))
    cv2.imwrite('tempflatfig.png', resized)

    zoomed_fig.savefig('tempfig.png')

    img = cv2.imread('tempfig.png')
    resized = cv2.resize(img, (600,300))
    cv2.imwrite('tempfig.png', resized)

    print('done')



import matplotlib
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import PySimpleGUI as sg

# Define the window layout
whole_graph = sg.Graph(canvas_size=(600,300), graph_bottom_left=(0, 0), graph_top_right=(600,300), enable_events=True, drag_submits=True, key='WholeGraph')
cmaps = ['oleron', 'vic', 'gray', 'acton', 'oslo', 'roma']

layout = [
    [sg.Text("Whole Moon Surface")],
    [sg.Text("Draw a rectangle over a region to zoom in")],
    [whole_graph, sg.Graph(canvas_size=(600,300), graph_bottom_left=(0, 0), graph_top_right=(600,300), enable_events=True, drag_submits=True, key='ZoomedGraph')],
    [sg.Text('Change the height of the 3D view',font=('Arial Bold', 12)), 
        sg.Button("Higher"), sg.Button('Lower')
    ],
    [sg.Text("Left: Zoomed In Area"), sg.Text("Right: Flat Version of Zoomed In Area")],
    [
     sg.Graph(canvas_size=(600,300), graph_bottom_left=(0, 0), graph_top_right=(600,300), enable_events=True, drag_submits=True, key='ZoomedFlatGraph')
     ],
    [sg.Text("Click on the flat version of the zoomed in area to sample point heights.")],
    [sg.Text('Longitude: '+"{:.3f}".format(selected_long)+' Latitude: '+"{:.3f}".format(selected_lat)+' Height: '+"{:.3f}".format(float(selected_point_height))+'m',font=('Arial Bold', 20), key='ChosenPointHeight', justification='left')]
]

# Create the form and show it without the plot
window = sg.Window(
    "General Public Moon Visualisation",
    layout,
    location=(0, 0),
    finalize=True,
    element_justification="center",
    font="Helvetica 18",
)

window['WholeGraph'].draw_image(filename='./wholefig.png', location=(0,300))
window['ZoomedGraph'].draw_image(filename='./tempfig.png', location=(0,300))
window['ZoomedFlatGraph'].draw_image(filename='./tempflatfig.png', location=(0,300))
minlong, maxlong, minlat, maxlat = region.split('/')

xmin = (float(minlong)/0.6) + 300
xmax = (float(maxlong)/0.6) + 300
ymin = (float(minlat)/0.6) + 150 
ymax = (float(maxlat)/0.6) + 150 

starting_rect = whole_graph.draw_rectangle((float(xmin), float(ymin)), (float(xmax), float(ymax)), line_color='red')

dragging = False
start_point = end_point = prior_rect = None

while True:
   event, values = window.read()
   if event == 'Higher':
      azimuth, elevation = perspective[0], perspective[1]
      azimuth = int(azimuth)
      elevation = int(elevation) + 10
      
      perspective = [azimuth, elevation]
      change_view()

      window['ZoomedGraph'].draw_image(filename='./tempfig.png', location=(0,300))

      window.refresh()

   if event == 'Lower':
      azimuth, elevation = perspective[0], perspective[1]
      azimuth = int(azimuth)
      elevation = int(elevation) - 10

      if elevation <= 15:
        elevation = 15
      
      perspective = [azimuth, elevation]
      change_view()

      window['ZoomedGraph'].draw_image(filename='./tempfig.png', location=(0,300))

      window.refresh()

   elif event == "WholeGraph":  # if there's a "Graph" event, then it's a mouse
        print('click: ', values[event])
        whole_graph.delete_figure(starting_rect)
        x, y = values[event]
        if not dragging:
            start_point = (x, y)
            dragging = True
        else:
            end_point = (x, y)
        if prior_rect:
            whole_graph.delete_figure(prior_rect)
        if None not in (start_point, end_point):
            prior_rect = whole_graph.draw_rectangle(start_point, end_point, line_color='red')

   elif event == 'ZoomedFlatGraph':
       x, y = values[event] 

       window['ZoomedFlatGraph'].erase()
       window['ZoomedFlatGraph'].draw_image(filename='./tempflatfig.png', location=(0,300))
       window['ZoomedFlatGraph'].draw_point((x,y), size=15, color='red')

       xmin, xmax, ymin, ymax = region.split('/')

       xmin, xmax, ymin, ymax = float(xmin), float(xmax), float(ymin), float(ymax) 

       actual_width = xmax - xmin
       actual_height = ymax - ymin

       long = (x / 600) * actual_width + xmin 
       lat = (y / 300) * actual_height + ymin 

       def find_nearest(array, value):
          array = np.asarray(array)
          idx = (np.abs(array - value)).argmin()
          return array[idx]
       
       selected_long = long 
       selected_lat = lat

       selected_point_height = new_grid.sel(x = selected_long, y = selected_lat, method = 'nearest').data

       window['ChosenPointHeight'].update('Longitude: '+"{:.3f}".format(selected_long)+' Latitude: '+"{:.3f}".format(selected_lat)+' Height: '+"{:.3f}".format(float(selected_point_height))+'m')


   elif event == 'WholeGraph+UP':  
       if start_point is not None and end_point is not None:
        start_lat = start_point[1] - 150
        start_long = start_point[0] - 300
            
        start_lat = 0.6 * start_lat
        start_long = 0.6 * start_long 

        end_lat = end_point[1] - 150
        end_long = end_point[0] - 300
            
        end_lat = 0.6 * end_lat
        end_long = 0.6 * end_long 

            # xmin/xmax/ymin/ymax
        lats = [start_lat, end_lat]
        longs = [start_long, end_long]
        xmin = min(longs)
        xmax = max(longs)
        ymin = min(lats)
        ymax = max(lats)

        minlong, maxlong, minlat, maxlat = float(xmin), float(xmax), float(ymin), float(ymax)

        region = str(xmin)+'/'+str(xmax)+'/'+str(ymin)+'/'+str(ymax)

        change_zoomed_region(region)
        window['ZoomedGraph'].draw_image(filename='./tempfig.png', location=(0,300))
        window['ZoomedFlatGraph'].draw_image(filename='./tempflatfig.png', location=(0,300))

        window.refresh()
        start_point, end_point = None, None 
        dragging = False

  
   if event == sg.WIN_CLOSED or event == 'Exit':
      break
window.close()