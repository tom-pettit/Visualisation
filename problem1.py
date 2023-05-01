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

# import the grid as an xarray DataArray and convert it to longitude and latitude
new_grid = xr.DataArray(
    data=dataarray,
    coords=dict(
        y=np.linspace(start=-90, stop=90, num=height),
        x=np.linspace(start=-180, stop=180, num=width),
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

# for if the scientist wants to view a specific point
selected_long = 0
selected_lat = 0
selected_point_height = new_grid[selected_long, selected_lat].data

# specify the height at which contour lines are drawn
contour_height = 1

# specify the default colour map
cmap = 'gray'

# create a grid image of the whole displacement map
whole_fig = pygmt.Figure()
whole_fig.grdimage(grid=new_grid, cmap=cmap, projection="Q12C")
whole_map_colorbar = whole_fig.colorbar(frame=["a1", "x+lElevation", "y+lkm"], position='JT')

# save the figure to a file
whole_fig.savefig('wholefig.png')

# create a grid view of the chosen zoomed in region
zoomed_fig = pygmt.Figure()
zoomed_fig.grdview(
    grid=new_grid,
    perspective=perspective,
    region=region,
    # sets the x- and y-axis labels, and annotates the west, south, and east axes
    frame=["xa", "ya", "WSnE"],
    # use the cylindrical equidistant projection on a figure of size 12cm 
    projection="Q12C",
    # sets the height of the three-dimensional relief at 2.5 centimeters
    zsize="2.5c",
    # apply the color map
    surftype="s",
    cmap=cmap
)

# create a colourmap for the zoomed in 3D perspective figure
zoomed_fig.colorbar(frame=["a1", "x+lElevation", "y+lkm"])
zoomed_fig.savefig('tempfig.png')

# create a figure for the flat version of the zoomed in region 
zoomed_flat_fig = pygmt.Figure()
zoomed_flat_fig.grdimage(grid=new_grid, cmap=cmap, projection="Q12C", region=region)
# add contour lines to the figure
zoomed_flat_fig.grdcontour(grid=new_grid, interval=1, annotation=1)

zoomed_flat_fig.savefig('tempflatfig.png')
import cv2 

# resize the images so that they can fit on the GUI (as the resize parameter from PyGMT does not work)
img = cv2.imread('wholefig.png')
resized = cv2.resize(img, (600,385))
cv2.imwrite('wholefig.png', resized)

img = cv2.imread('tempfig.png')
resized = cv2.resize(img, (600,300), interpolation=cv2.INTER_AREA)
cv2.imwrite('tempfig.png', resized)

img = cv2.imread('tempflatfig.png')
resized = cv2.resize(img, (600,300), interpolation=cv2.INTER_AREA)
cv2.imwrite('tempflatfig.png', resized)



# function to change the perspective of the 3D perspective image
def change_perspective(perspective):
    import os
    os.remove('tempfig.png')
    # extract the new azimuth and elevation values from the inputted perspective
    azimuth, elevation = perspective[1:-1].split(',')
    azimuth = int(azimuth)
    elevation = int(elevation)
    zoomed_fig = pygmt.Figure()
    zoomed_fig.grdview(
        grid=new_grid,
        # sets the new perspective
        perspective=[azimuth, elevation],
        region=region,
        frame=["xa", "ya", "WSnE"],
        projection="Q12C",
        zsize="2.5c",
        surftype="s",
        cmap=cmap
    )
    
    # the color bar
    zoomed_fig.colorbar(frame=["a1", "x+lElevation", "y+lkm"])

    zoomed_fig.savefig('tempfig.png')

    # resize the image for the GUI
    img = cv2.imread('tempfig.png')
    resized = cv2.resize(img, (600,300), interpolation=cv2.INTER_AREA)
    cv2.imwrite('tempfig.png', resized)

# function to change the region the scientist would like to zoom into
def change_zoomed_region(region):
    import os
    os.remove('tempfig.png')
    
    zoomed_fig = pygmt.Figure()
    zoomed_fig.grdview(
        grid=new_grid,
        perspective=perspective,
        # set the region as the new chosen region
        region=region,
        frame=["xa", "ya", "WSnE"],
        projection="Q12C",
        zsize="2.5c",
        surftype="s",
        cmap=cmap
    )

    # the colorbar
    zoomed_fig.colorbar(frame=["a1", "x+lElevation", "y+lkm"])

    # re-create the flat version of the zoomed in region
    zoomed_flat_fig = pygmt.Figure()
    zoomed_flat_fig.grdimage(grid=new_grid, cmap=cmap, projection="Q12C", region=region)
    # add the contour lines
    zoomed_flat_fig.grdcontour(grid=new_grid, interval=1, annotation=1)

    # output the flat version of the 3D perspective figure to a file and resize for the GUI
    zoomed_flat_fig.savefig('tempflatfig.png')
    img = cv2.imread('tempflatfig.png')
    resized = cv2.resize(img, (600,300), interpolation=cv2.INTER_AREA)
    cv2.imwrite('tempflatfig.png', resized)

    # output the 3D perspective figure to a file and resize for the GUI
    zoomed_fig.savefig('tempfig.png')
    img = cv2.imread('tempfig.png')
    resized = cv2.resize(img, (600,300), interpolation=cv2.INTER_AREA)
    cv2.imwrite('tempfig.png', resized)

# function to change the height at which contour lines are drawn on the flat version of the region
def change_contours(height_choice):
    import os
    os.remove('tempflatfig.png')

    # re-create the flat version of the 3D perspective figure
    zoomed_flat_fig = pygmt.Figure()
    zoomed_flat_fig.grdimage(grid=new_grid, cmap=cmap, projection="Q12C", region=region)
    # add the new contour lines with the inputted height
    zoomed_flat_fig.grdcontour(grid=new_grid, interval=height_choice, annotation=1)

    # save the figure to a file and resize for the GUI
    zoomed_flat_fig.savefig('tempflatfig.png')
    img = cv2.imread('tempflatfig.png')
    resized = cv2.resize(img, (600,300), interpolation=cv2.INTER_AREA)
    cv2.imwrite('tempflatfig.png', resized)

# function to change the colour map for all the figures on the GUI
def change_colourmap(cmap):
    # create the figure for the whole surface of the moon
    whole_fig = pygmt.Figure()
    whole_fig.grdimage(grid=new_grid, cmap=cmap, projection="Q12C")
    # add the colorbar for this new figure
    whole_map_colorbar = whole_fig.colorbar(frame=["a1", "x+lElevation", "y+lkm"], position='JT')

    # save the figure to a file
    whole_fig.savefig('wholefig.png')

    # re-create the 3D perspective figure
    zoomed_fig = pygmt.Figure()
    zoomed_fig.grdview(
        grid=new_grid,
        perspective=perspective,
        region=region,
        frame=["xa", "ya", "WSnE"],
        projection="Q12C",
        zsize="4.5c",
        surftype="s",
        # sets the colour map as the new chosen colour map 
        cmap=cmap
    )
    # add the colorbar for the new figure
    zoomed_fig.colorbar(frame=["a1", "x+lElevation", "y+lkm"])
    zoomed_fig.savefig('tempfig.png')

    # re-create the figure for the flat version of the 3D perspective figure
    zoomed_flat_fig = pygmt.Figure()
    zoomed_flat_fig.grdimage(grid=new_grid, cmap=cmap, projection="Q12C", region=region)
    zoomed_flat_fig.grdcontour(grid=new_grid, interval=1, annotation=1)
    zoomed_flat_fig.savefig('tempflatfig.png')

    # resize all the newly created figures for the GUI
    img = cv2.imread('wholefig.png')
    resized = cv2.resize(img, (600,385))
    cv2.imwrite('wholefig.png', resized)

    img = cv2.imread('tempflatfig.png')
    resized = cv2.resize(img, (600,300), interpolation=cv2.INTER_AREA)
    cv2.imwrite('tempflatfig.png', resized)

    img = cv2.imread('tempfig.png')
    resized = cv2.resize(img, (600,300), interpolation=cv2.INTER_AREA)
    cv2.imwrite('tempfig.png', resized)

# use PySimpleGUI for the GUI
import PySimpleGUI as sg

# Define the window layout

# the object for the whole moon surface figure
whole_graph = sg.Graph(canvas_size=(600,385), graph_bottom_left=(0, 0), graph_top_right=(600,385), enable_events=True, drag_submits=True, key='WholeGraph')

# the options for the color map
cmaps = ['oleron', 'vic', 'gray', 'acton', 'magma', 'roma']

# the dropdown list object allowing a scientist to change the colour map of the figures
lst = sg.Combo(cmaps, font=('Arial Bold', 14),  default_value='gray', expand_x=True, enable_events=True,  readonly=False, key='CmapDropdown')

# the layout array for the GUI
layout = [
    [sg.Text("Whole Moon Surface")],
    [sg.Text("Draw a rectangle over a region to zoom in")],
    [whole_graph, lst, sg.Button('Change Colourmap')],
    [sg.Text("Viewing Moon Surface at Region: "+"Long: "+str(minlong)+" to "+str(maxlong)+" Lat: "+str(minlat)+" to "+str(maxlat), key='region_title'),
     sg.Text('At Azimuth: '+str(perspective[0])+' and Elevation: '+str(perspective[1]), key='perspective_title'),
     sg.Text('Contour Height Display Interval: '+str(contour_height), key='contour_height')],
    [sg.Text("Left: Zoomed In Area"), sg.Text("Right: Flat Version of Zoomed In Area")],
    [sg.Graph(canvas_size=(600,300), graph_bottom_left=(0, 0), graph_top_right=(600,300), enable_events=True, drag_submits=True, key='ZoomedGraph'),
     sg.Graph(canvas_size=(600,300), graph_bottom_left=(0, 0), graph_top_right=(600,300), enable_events=True, drag_submits=True, key='ZoomedFlatGraph')
     ],
    [sg.Text("Click on the flat version of the zoomed in area to sample point heights.")],
    [sg.Text('Longitude: '+"{:.3f}".format(selected_long)+' Latitude: '+"{:.3f}".format(selected_lat)+' Height: '+"{:.3f}".format(float(selected_point_height))+'km',font=('Arial Bold', 20), key='ChosenPointHeight', justification='left')],

    [sg.Text('Enter a View Azimuth and View Elevation in the form: [View Azimuth, View Elevation]',font=('Arial Bold', 12)), 
        sg.Input('', enable_events=True, key='Perspective', font=('Arial Bold', 20), expand_x=True, justification='center'),
        sg.Button("Change Perspective")
    ],
    [sg.Text('Enter a Height at which Isocontours are Generated (Default: 1): ',font=('Arial Bold', 12)), 
        sg.Input('', enable_events=True, key='Iscontour_height', font=('Arial Bold', 20), expand_x=True, justification='center'),
        sg.Button("Change Isocontours")
    ]
]

# Create the window object for the GUI
window = sg.Window(
    # title
    "Scientific Moon Visualisation",
    # set the layout as the layout specified above
    layout,
    location=(0, 0),
    finalize=True,
    element_justification="center",
    font="Helvetica 18",
)

# draw the PyGMT figures on the GUI
window['WholeGraph'].draw_image(filename='./wholefig.png', location=(0,385))
window['ZoomedGraph'].draw_image(filename='./tempfig.png', location=(0,300))
window['ZoomedFlatGraph'].draw_image(filename='./tempflatfig.png', location=(0,300))

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
   event, values = window.read()
   # if the user presses the change perspective button
   if event == 'Change Perspective':
      new_perspective = values['Perspective']
      # call the change perspective function
      change_perspective(new_perspective)
      # re-draw the 3D perspective figure
      window['ZoomedGraph'].draw_image(filename='./tempfig.png', location=(0,300))
      azimuth, elevation = new_perspective[1:-1].split(',')
      azimuth = int(azimuth)
      elevation = int(elevation)
      perspective = [azimuth, elevation]
      # update the information on the GUI representing the perspective of the 3D perspective figure
      window['perspective_title'].update('At Azimuth: '+str(perspective[0])+' and Elevation: '+str(perspective[1]))
      # refresh the window
      window.refresh()

   # if the user presses the change isocontours button
   if event == 'Change Isocontours':
      # read in the newly chosen isocontour height
      new_height = float(values['Iscontour_height'])
      contour_height = new_height
      # call the change contours function
      change_contours(new_height)
      # re-draw the new flat zoomed in figure with the new isocontours
      window['ZoomedFlatGraph'].draw_image(filename='./tempflatfig.png', location=(0,300))
      # update the text representing the current isocontour settings
      window['contour_height'].update('Contour Height Display Interval: '+str(contour_height))
      # refresh the window
      window.refresh()

   # if the user presses the change colourmap button
   if event == 'Change Colourmap':
      # read the chosen colormap from the dropdown menu
      cmap = values['CmapDropdown']
      # call the change colormap function
      change_colourmap(cmap)

      # re-draw all 3 figures on the GUI with the new colormap
      window['WholeGraph'].draw_image(filename='./wholefig.png', location=(0,385))
      window['ZoomedGraph'].draw_image(filename='./tempfig.png', location=(0,300))
      window['ZoomedFlatGraph'].draw_image(filename='./tempflatfig.png', location=(0,300))

      # calculate where the rectangle should be drawn
      minlong, maxlong, minlat, maxlat = region.split('/')

      xmin = (float(minlong)/0.6) + 300
      xmax = (float(maxlong)/0.6) + 300
      ymin = (float(minlat)/0.6) + 150
      ymax = (float(maxlat)/0.6) + 150

      # re-draw the red rectangle
      starting_rect = whole_graph.draw_rectangle((float(xmin), float(ymin)), (float(xmax), float(ymax)), line_color='red')

      # refresh the window
      window.refresh()

   # if the user interacts with the whole surface of the moon figure
   elif event == "WholeGraph":  
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

   # if the user clicks on the flat version of the zoomed in region
   elif event == 'ZoomedFlatGraph':
       # read the mouse position
       x, y = values[event] 

       # erase the figure
       window['ZoomedFlatGraph'].erase()
       # re-draw the figure
       window['ZoomedFlatGraph'].draw_image(filename='./tempflatfig.png', location=(0,300))
       # add a red dot at the location of the scientist's click
       window['ZoomedFlatGraph'].draw_point((x,y), size=15, color='red')

       # calculate the longitude and latitude values from the coordinates of the mouse click
       xmin, xmax, ymin, ymax = region.split('/')

       xmin, xmax, ymin, ymax = float(xmin), float(xmax), float(ymin), float(ymax) 

       actual_width = xmax - xmin
       actual_height = ymax - ymin

       long = (x / 600) * actual_width + xmin 
       lat = (y / 300) * actual_height + ymin 
       
       selected_long = long 
       selected_lat = lat

       # collect the elevation of the chosen point, using the nearest point in the dataarray
       # due to the sampling rate of the grid, the scientist may not pick an exact datapoint, so select the nearest
       selected_point_height = new_grid.sel(x = selected_long, y = selected_lat, method = 'nearest').data

       # udpate the text showing the scientist the point that they clicked and the height of the surface at that point
       window['ChosenPointHeight'].update('Longitude: '+"{:.3f}".format(selected_long)+' Latitude: '+"{:.3f}".format(selected_lat)+' Height: '+"{:.3f}".format(float(selected_point_height))+'km')

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

        # re-draw the 3D perspective figure and the flat version of it
        window['ZoomedGraph'].draw_image(filename='./tempfig.png', location=(0,300))
        window['ZoomedFlatGraph'].draw_image(filename='./tempflatfig.png', location=(0,300))

        # update the text showing the scientist the current settings of their zoomed in region
        window['region_title'].update("Viewing Moon Surface at Region: "+"Long: "+"{:.3f}".format(minlong)+" to "+"{:.3f}".format(maxlong)+" Lat: "+"{:.3f}".format(minlat)+" to "+"{:.3f}".format(maxlat))

        # refresh the window
        window.refresh()

        # reset the variables to help with drawing the rectangle
        start_point, end_point = None, None 
        dragging = False

   # if the scientist wants to close the window, exit the loop
   if event == sg.WIN_CLOSED or event == 'Exit':
      break
   
# close the window if the loop is broken
window.close()