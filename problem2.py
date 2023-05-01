import pygmt
import matplotlib.pyplot as plt
import rioxarray

# load in the displacement dataset
dataarray =  rioxarray.open_rasterio('./dataset/displacement_dataset/ldem_4.tif')
tif_shape = dataarray.shape
dataarray = dataarray.isel(band=0)
import xarray as xr
import numpy as np 

height, width = dataarray.shape

xprecision = 360 / width 
yprecision = 180 / height

# create an xarray DataArray for the displacement dataset
new_grid = xr.DataArray(
    data=dataarray,
    coords=dict(
        y=np.linspace(start=-90, stop=90, num=height),
        x=np.linspace(start=-180, stop=180, num=width),
    ),
    dims=("y", "x"),
)

# new_grid = pygmt.grdsample(
#     grid=new_grid, translate=True, spacing=[0.25, 0.25]
# )

# load in the colour dataset
color_dataarray =  rioxarray.open_rasterio('./dataset/color_dataset/lroc_color_poles_4k.tif')
color_tif_shape = color_dataarray.shape

color_dataarray = color_dataarray.isel(band=0)
import xarray as xr
import numpy as np 

color_height, color_width = color_dataarray.shape
# print(color_dataarray)

# create an xarray DataArray for the color dataset
color_grid = xr.DataArray(
    data=color_dataarray,
    coords=dict(
        y=np.linspace(start=0, stop=color_height, num=color_height),
        x=np.linspace(start=0, stop=color_width, num=color_width),
    ),
    dims=("y", "x"),
)

# specify some parameters for the visualisation
# specify the region of the zoomed in area
region = "90/110/5/25"
# convert the region string into longitude and latitude values in order to draw a rectangle
minlong, maxlong, minlat, maxlat = region.split('/')
minlong, maxlong, minlat, maxlat = float(minlong), float(maxlong), float(minlat), float(maxlat)
perspective = [150, 30]

# for if the user wants to view a specific point
selected_long = 0
selected_lat = 0
selected_point_height = new_grid[selected_long, selected_lat].data

# specify the default colour map
cmap = 'gray'

# create a figure for the whole moon surface using the colour map 
whole_fig = pygmt.Figure()
color_grid = pygmt.grdproject(grid=color_grid, projection="x1:400")
whole_fig.grdimage(grid=color_grid, cmap='gray')

# save the whole moon surface figure to a file
whole_fig.savefig('wholefig.png')

# create a figure for the 3D perspective image
zoomed_fig = pygmt.Figure()
zoomed_fig.grdview(
    grid=new_grid,
    perspective=perspective,
    region=region,
    # sets the x- and y-axis labels, and annotates the west, south, and east axes
    frame=["xa", "ya", "WSnE"],
    # use the cylindrical equidistant projection on a figure of size 12cm
    projection="Q12C",
    # sets the height of the three-dimensional relief at 4.5 centimeters
    zsize="4.5c",
    # apply the color map 
    surftype="s",
    cmap=cmap
)

# create a colourmap for the zoomed in 3D perspective figure
zoomed_fig.colorbar(frame=["a1", "x+lElevation", "y+lkm"])
zoomed_fig.savefig('tempfig.png')
import cv2 

# resize the figures to fit in the GUI
img = cv2.imread('wholefig.png')
resized = cv2.resize(img, (600,300))
cv2.imwrite('wholefig.png', resized)

img = cv2.imread('tempfig.png')
resized = cv2.resize(img, (550,300), interpolation=cv2.INTER_AREA)
cv2.imwrite('tempfig.png', resized)

# create the histogram figure for showing the heights of possible landing sites
fig = pygmt.Figure()
fig.histogram(
    data=dataarray.data,
    # define the frame, add a title, and set the background color to
    # "lightgray". Add labels to the x-axis and y-axis
    frame=["WSne+tHistogram+glightgray", "x+lElevation (km)", "y+lCounts"],
    # bins of size 1 (for each metre of elevation)
    series=1,
    # add numbers to the bars
    annotate=True,
    # specify the Artemis 3 colours for the bar fills
    fill="#060b27",
    # use the pen parameter to draw the outlines with a width of 1 point
    pen="1p",
    # histogram type is counts
    histtype=0,
)

# save the histogram figure to a file
fig.savefig('histogram.png')

# resize the histogram to fit on the GUI
img = cv2.imread('histogram.png')
resized = cv2.resize(img, (400,300), interpolation=cv2.INTER_AREA)
cv2.imwrite('histogram.png', resized)

# resize the Artemis 3 mission logo to fit on the top of the GUI
img = cv2.imread('logo.png')
resized = cv2.resize(img, (300,100))
cv2.imwrite('logo_resized.png', resized)

# function to change the view of the 3D perspective image
def change_view():
    import os
    os.remove('tempfig.png')
    # get the azimuth and elevation 
    azimuth, elevation = perspective[0], perspective[1]
    azimuth = int(azimuth)
    elevation = int(elevation)

    # create the 3D perspective figure
    zoomed_fig = pygmt.Figure()
    zoomed_fig.grdview(
        grid=new_grid,
        # set the perspective at the new perspective 
        perspective=[azimuth, elevation],
        region=region,
        frame=["xa", "ya", "WSnE"],
        projection="Q12C",
        zsize="4.5c",
        surftype="s",
        cmap=cmap
    )
    
    # add the color bar 
    zoomed_fig.colorbar(frame=["a1", "x+lElevation", "y+lkm"])

    # save the figure to a file
    zoomed_fig.savefig('tempfig.png')

    # resize the figure to fit on the GUI
    img = cv2.imread('tempfig.png')
    resized = cv2.resize(img, (550,300), interpolation=cv2.INTER_AREA)
    cv2.imwrite('tempfig.png', resized)

# function to change the region of the 3D perspective figure
def change_zoomed_region(region):
    import os
    os.remove('tempfig.png')
    
    # create the figure
    zoomed_fig = pygmt.Figure()
    zoomed_fig.grdview(
        grid=new_grid,
        perspective=perspective,
        # set the region as the new region
        region=region,
        frame=["xa", "ya", "WSnE"],
        projection="Q12C",
        zsize="4.5c",
        surftype="s",
        cmap=cmap
    )

    # add the colour bar
    zoomed_fig.colorbar(frame=["a1", "x+lElevation", "y+lkm"])

    # save the figure to a file
    zoomed_fig.savefig('tempfig.png')

    # resize the figure to fit on the GUI
    img = cv2.imread('tempfig.png')
    resized = cv2.resize(img, (550,300), interpolation=cv2.INTER_AREA)
    cv2.imwrite('tempfig.png', resized)

# use PySimpleGUI to create the GUI
import PySimpleGUI as sg

# Define the window layout

# the object for the whole moon surface figure
whole_graph = sg.Graph(canvas_size=(600,300), graph_bottom_left=(0, 0), graph_top_right=(600,300), enable_events=True, drag_submits=True, key='WholeGraph')

# the options for the color map
cmaps = ['oleron', 'vic', 'gray', 'acton', 'oslo', 'roma']

# the layout array for the GUI
layout = [
    [sg.Image('./rotating_moon2.gif', key = "MoonGIF", background_color='#060b27'), sg.Image('./logo_resized.png', background_color='#060b27')],
    [sg.Text("Combatting the complex terrain, to take the next wave of humans to the moon.", font=('Avenir', 25), background_color='#060b27')],
    [sg.Text("The picture below shows the surface of our moon. Its uneven terrain is tricky to land a rocket. Have a look for yourself, by drawing a rectangle over a region to zoom in!", font=('Avenir', 15), background_color='#060b27')],
    [sg.Text("On the far right is the heights of all the possible landing sites. Can you find them?", font=('Avenir', 15), background_color='#060b27')],
    [whole_graph, sg.Graph(canvas_size=(550,300), graph_bottom_left=(0, 0), graph_top_right=(550,300), enable_events=True, drag_submits=True, key='ZoomedGraph'), sg.Image('./histogram.png', background_color='#060b27')],
    [sg.Text('Change the height of the 3D view',font=('Avenir', 12), background_color='#060b27'), 
        sg.Button("Higher"), sg.Button('Lower'), sg.Button('Left'), sg.Button('Right')
    ],
    [sg.Text("Why are we doing this?", font=('Avenir', 30), background_color='#060b27')],
    [sg.Text("1. We’re building on more than 50 years of exploration experience to reignite America’s passion for discovery.", font=('Avenir', 18), background_color='#060b27')],
    [sg.Text("2. Artemis missions enable a growing lunar economy by fueling new industries, supporting job growth, and furthering the demand for a skilled workforce.", font=('Avenir', 18), background_color='#060b27')],
    [sg.Text("3. Inspiring the next generation", font=('Avenir', 18), background_color='#060b27')]
    
]

# Create the window object for the GUI
window = sg.Window(
    # title
    "General Public Moon Visualisation",
    # set the layout as the layout specified above
    layout,
    location=(0, 0),
    finalize=True,
    element_justification="center",
    font="Avenir 18",
    background_color='#060b27'
)

# draw the PyGMT figures on the GUI
window['WholeGraph'].draw_image(filename='./wholefig.png', location=(0,300))
window['ZoomedGraph'].draw_image(filename='./tempfig.png', location=(0,300))

# calculate the rectangle to draw on the whole moon surface figure representing the area that the 3D perspective image is at by default
minlong, maxlong, minlat, maxlat = region.split('/')

xmin = (float(minlong)/0.6) + 300
xmax = (float(maxlong)/0.6) + 300
ymin = (float(minlat)/0.6) + 150 
ymax = (float(maxlat)/0.6) + 150 

# draw the starting rectangle for the originally selected region of the moon shown in the 3D perspective figure
starting_rect = whole_graph.draw_rectangle((float(xmin), float(ymin)), (float(xmax), float(ymax)), line_color='red')

# variables to help draw the rectangle specified by the scientist
dragging = False
start_point = end_point = prior_rect = None

# create a loop that can read events by the scientist when they interract with the GUI
while True:
    # read the events taking place on the window
   event, values = window.read(timeout=50)
   # if the user has pressed the higher button
   if event == 'Higher':
      azimuth, elevation = perspective[0], perspective[1]
      azimuth = int(azimuth)

      # increase the elevation of the 3D perspective figure by 10 metres
      elevation = int(elevation) + 10
      
      # re-create the perspective variable
      perspective = [azimuth, elevation]

      # call the change view function
      change_view()

      # re-draw the 3D perspective figure on the GUI
      window['ZoomedGraph'].draw_image(filename='./tempfig.png', location=(0,300))

      # refresh the window
      window.refresh()

   # if thre user has pressed teh lower button
   if event == 'Lower':
      azimuth, elevation = perspective[0], perspective[1]
      azimuth = int(azimuth)

      # decrease the elevation of the 3D perspective figure by 10 metres
      elevation = int(elevation) - 10

      # ensure the elevation of the 3D perspective image does not go lower than 15 metres
      if elevation <= 15:
        elevation = 15
      
      # re-create the perspective variable
      perspective = [azimuth, elevation]

      # call the change view function 
      change_view()

      # re-draw the 3D perspective figure on the GUI
      window['ZoomedGraph'].draw_image(filename='./tempfig.png', location=(0,300))

      # refresh the window
      window.refresh()

   # if the user presses the left button
   if event == 'Left':
      azimuth, elevation = perspective[0], perspective[1]
      # rotate the view azimuth by 20 degrees in the leftwards direction
      azimuth = int(azimuth) + 20
      elevation = int(elevation) 
      
      # re-create the perspective variable
      perspective = [azimuth, elevation]

      # call the change view function
      change_view()

      # re-draw the 3D perspective image on the GUI
      window['ZoomedGraph'].draw_image(filename='./tempfig.png', location=(0,300))

      # refresh the window
      window.refresh()

   # if the user presses the right button
   if event == 'Right':
      azimuth, elevation = perspective[0], perspective[1]
      # rotate the view azimuth by 20 degrees in the rightwards direction
      azimuth = int(azimuth) - 20
      elevation = int(elevation) 
      
      # re-create the perspective variable
      perspective = [azimuth, elevation]

      # call the change view function
      change_view()

      # re-draw the 3D perspective figure on the GUI
      window['ZoomedGraph'].draw_image(filename='./tempfig.png', location=(0,300))

      # refresh the window
      window.refresh()

   # if the user interacts with the whole surface of the moon figure
   elif event == "WholeGraph":  # if there's a "Graph" event, then it's a mouse
        # delete the current rectangle drawn on the figure
        whole_graph.delete_figure(starting_rect)
        # get the coordinates of the mouse clicking on the figure
        x, y = values[event]
        # if the user has just clicked on the figure for the first time
        if not dragging:
            start_point = (x, y)
            dragging = True
        else:
            # if the user has finished drawing their rectangle
            end_point = (x, y)

        # delete any previous rectangle 
        if prior_rect:
            whole_graph.delete_figure(prior_rect)

        # if both the start and end point of the rectangle have been specified
        if None not in (start_point, end_point):
            # draw the rectangle
            prior_rect = whole_graph.draw_rectangle(start_point, end_point, line_color='red')

   # if the scientist has stopped dragging the rectangular box for the zoomed in region 
   elif event == 'WholeGraph+UP':  
        # if the start and end points of the rectangle have been received
       if start_point is not None and end_point is not None:
        # convert the coordinates into longitude and latitude values
        start_lat = start_point[1] - 150
        start_long = start_point[0] - 300
            
        start_lat = 0.6 * start_lat
        start_long = 0.6 * start_long 

        end_lat = end_point[1] - 150
        end_long = end_point[0] - 300
            
        end_lat = 0.6 * end_lat
        end_long = 0.6 * end_long 

        # update the region variable in the correct format
        # ie the xmin, xmax, ymin, ymax need to be in the correct order
        lats = [start_lat, end_lat]
        longs = [start_long, end_long]
        xmin = min(longs)
        xmax = max(longs)
        ymin = min(lats)
        ymax = max(lats)

        minlong, maxlong, minlat, maxlat = float(xmin), float(xmax), float(ymin), float(ymax)
        region = str(xmin)+'/'+str(xmax)+'/'+str(ymin)+'/'+str(ymax)

        # call the changed zoomed region function
        change_zoomed_region(region)

        # re-draw the 3D perspective figure 
        window['ZoomedGraph'].draw_image(filename='./tempfig.png', location=(0,300))

        # refresh the window
        window.refresh()

        # reset the variables to help with drawing the rectangle
        start_point, end_point = None, None 
        dragging = False

   # update the GIF of the moon so that it rotates
   window["MoonGIF"].UpdateAnimation("./rotating_moon2.gif",time_between_frames=50)

   # if the scientist wants to close the window, exit the loop
   if event == sg.WIN_CLOSED or event == 'Exit':
      break
   
# close the window if the loop is broken
window.close()