#region imports
from scipy.integrate import odeint
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import numpy as np
import math
from PyQt5 import QtWidgets as qtw
from PyQt5 import QtCore as qtc
from PyQt5 import QtGui as qtg

#these imports are necessary for drawing a matplot lib graph on my GUI
#no simple widget for this exists in QT Designer, so I have to add the widget in code.
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure
#endregion

#region class definitions
#region specialized graphic items
class MassBlock(qtw.QGraphicsItem):
    def __init__(self, CenterX, CenterY, width=30, height=10, parent=None, pen=None, brush=None, name='CarBody', mass=10):
        super().__init__(parent)
        self.x = CenterX
        self.y = CenterY
        self.pen = pen
        self.brush = brush
        self.width = width
        self.height = height
        self.top = self.y - self.height/2
        self.left = self.x - self.width/2
        self.rect = qtc.QRectF(self.left, self.top, self.width, self.height)
        self.name = name
        self.mass = mass
        self.transformation = qtg.QTransform()
        stTT = self.name +"\nx={:0.3f}, y={:0.3f}\nmass = {:0.3f}".format(self.x, self.y, self.mass)
        self.setToolTip(stTT)

    def boundingRect(self):
        bounding_rect = self.transformation.mapRect(self.rect)
        return bounding_rect

    def paint(self, painter, option, widget=None):
        self.transformation.reset()
        if self.pen is not None:
            painter.setPen(self.pen)  # Red color pen
        if self.brush is not None:
            painter.setBrush(self.brush)
        self.top = -self.height/2
        self.left = -self.width/2
        self.rect=qtc.QRectF( self.left, self.top, self.width, self.height)
        painter.drawRect(self.rect)
        self.transformation.translate(self.x, self.y)
        self.setTransform(self.transformation)
        self.transformation.reset()
        # brPen=qtg.QPen()
        # brPen.setWidth(0)
        # painter.setPen(brPen)
        # painter.setBrush(qtc.Qt.NoBrush)
        # painter.drawRect(self.boundingRect())

class Wheel(qtw.QGraphicsItem):
    def __init__(self, CenterX, CenterY, radius=10, parent=None, pen=None, wheelBrush=None, massBrush=None, name='Wheel', mass=10):
        super().__init__(parent)
        self.x = CenterX
        self.y = CenterY
        self.pen = pen
        self.brush = wheelBrush
        self.radius = radius
        self.rect = qtc.QRectF(self.x - self.radius, self.y - self.radius, self.radius*2, self.radius*2)
        self.name = name
        self.mass = mass
        self.transformation = qtg.QTransform()
        stTT = self.name +"\nx={:0.3f}, y={:0.3f}\nmass = {:0.3f}".format(self.x, self.y, self.mass)
        self.setToolTip(stTT)
        self.massBlock = MassBlock(CenterX, CenterY, width=2*radius*0.85, height=radius/3, pen=pen, brush=massBrush, name="Wheel Mass", mass=mass)

    def boundingRect(self):
        bounding_rect = self.transformation.mapRect(self.rect)
        return bounding_rect
    def addToScene(self, scene):
        scene.addItem(self)
        scene.addItem(self.massBlock)

    def paint(self, painter, option, widget=None):
        self.transformation.reset()
        if self.pen is not None:
            painter.setPen(self.pen)  # Red color pen
        if self.brush is not None:
            painter.setBrush(self.brush)
        self.rect=qtc.QRectF(-self.radius, -self.radius, self.radius*2, self.radius*2)
        painter.drawEllipse(self.rect)
        self.transformation.translate(self.x, self.y)
        self.setTransform(self.transformation)
        self.transformation.reset()
        # brPen=qtg.QPen()
        # brPen.setWidth(0)
        # painter.setPen(brPen)
        # painter.setBrush(qtc.Qt.NoBrush)
        # painter.drawRect(self.boundingRect())

class Spring(qtw.QGraphicsItem):
    def __init__(self, x1, y1, x2, y2, parent=None, pen=None):
        super().__init__(parent)  # Initialize parent QGraphicsItem
        self.x1 = x1  # Starting x-coordinate
        self.y1 = y1  # Starting y-coordinate
        self.x2 = x2  # Ending x-coordinate
        self.y2 = y2  # Ending y-coordinate
        self.pen = pen if pen else qtg.QPen(qtc.Qt.black)  # Default pen if none provided
        self.updateEndpoints(x1, y1, x2, y2)  # Update geometry based on endpoints

    def boundingRect(self):
        # Create bounding rectangle expanded a little for drawing safety
        return qtc.QRectF(0, 0, self.width, self.height).adjusted(-10, -10, 10, 10)

    def paint(self, painter, option, widget=None):
        painter.setPen(self.pen)  # Set pen
        path = qtg.QPainterPath()  # Create a painter path
        path.moveTo(self.width / 2, 0)  # Start path at the middle-top
        n = 8  # Number of zigzag segments
        segment_height = self.height / (n + 2)  # Calculate vertical spacing for zigzag
        path.lineTo(self.width / 2, segment_height)  # Start vertical line
        # Create zigzag pattern by alternating sides
        for i in range(n):
            x = 0 if i % 2 == 0 else self.width  # Alternate x between left and right
            y = segment_height * (i + 2)  # Move down by segment height
            path.lineTo(x, y)
        path.lineTo(self.width / 2, self.height)  # Finish at center-bottom
        painter.drawPath(path)  # Draw the final spring shape

    def updateEndpoints(self, x1, y1, x2, y2):
        self.prepareGeometryChange()  # Prepare for changes in geometry
        self.setPos(x1, y1)  # Set position of spring
        self.x1 = x1  # Update internal x1
        self.y1 = y1  # Update internal y1
        self.x2 = x2  # Update internal x2
        self.y2 = y2  # Update internal y2
        self.width = 20  # Fixed width of spring drawing
        self.height = y2 - y1  # Set height based on start and end points
        self.update()  # Force item to redraw

# Class that defines a dashpot (shock absorber) graphical object
class DashpotItem(qtw.QGraphicsItem):
    def __init__(self, x1, y1, x2, y2, parent=None, pen=None):
        super().__init__(parent)  # Initialize parent QGraphicsItem
        self.x1 = x1  # Start x-coordinate
        self.y1 = y1  # Start y-coordinate
        self.x2 = x2  # End x-coordinate
        self.y2 = y2  # End y-coordinate
        self.pen = pen if pen else qtg.QPen(qtc.Qt.black)  # Default black pen if none provided
        self.updateEndpoints(x1, y1, x2, y2)  # Update geometry

    def boundingRect(self):
        # Create bounding rectangle for the dashpot
        return qtc.QRectF(0, 0, self.width, self.height).adjusted(-10, -10, 10, 10)

    def paint(self, painter, option, widget=None):
        painter.setPen(self.pen)  # Set pen for drawing
        painter.setBrush(qtg.QBrush(qtg.QColor(200, 200, 200, 128)))  # Light gray brush for filling the cylinder

        cylinder_height = self.height * 0.6  # Set 60% height for cylinder body
        cylinder_rect = qtc.QRectF(0, 0, self.width, cylinder_height)  # Define rectangle for cylinder
        painter.drawRect(cylinder_rect)  # Draw cylinder rectangle

        rod_width = self.width * 0.3  # Width of piston rod inside the cylinder
        rod_rect = qtc.QRectF(self.width / 2 - rod_width / 2, cylinder_height, rod_width, self.height - cylinder_height)  # Define rod rectangle
        painter.drawRect(rod_rect)  # Draw piston rod

        painter.drawLine(0, cylinder_height / 2, self.width, cylinder_height / 2)  # Draw a horizontal line across cylinder center

    def updateEndpoints(self, x1, y1, x2, y2):
        self.prepareGeometryChange()  # Prepare for geometry change
        self.setPos(x1, y1)  # Set new position
        self.x1 = x1  # Store new x1
        self.y1 = y1  # Store new y1
        self.x2 = x2  # Store new x2
        self.y2 = y2  # Store new y2
        self.width = 20  # Fixed width
        self.height = y2 - y1  # Update height based on endpoints
        self.update()  # Force redraw
#endregion

#endregion

#region MVC for quarter car model
#region MVC for quarter car model

# 🚗 Class that defines the "Model" — it stores all the car's properties and simulation results
class CarModel():
    """
    I re-wrote the quarter car model as an object-oriented program
    and used the MVC pattern. This class stores information about the car and
    results of the ODE calculation.
    """
    def __init__(self):
        # Simulation results from solving the ODEs (starts empty)
        self.results = []
        # Max simulation time in seconds
        self.tmax = 3.0
        # Create a time vector from 0 to tmax with 200 points
        self.t = np.linspace(0, self.tmax, 200)
        # Time required to climb the ramp (in seconds)
        self.tramp = 1.0
        # Ramp angle in radians (default small angle)
        self.angrad = 0.1
        # Ramp height in meters (6 inches divided by 12 inches/foot * 3.3 ft/m)
        self.ymag = 6.0 / (12 * 3.3)
        # Ramp angle in degrees (default 45° for a steep ramp)
        self.yangdeg = 45.0
        # Set results to None initially (until simulation is run)
        self.results = None

        # Set default physical parameters of the car model
        self.m1 = 450  # Mass of car body in kg
        self.m2 = 20  # Mass of wheel in kg
        self.c1 = 4500  # Damping coefficient of suspension (N·s/m)
        self.k1 = 15000  # Spring constant of suspension (N/m)
        self.k2 = 90000  # Spring constant of tire (N/m)
        self.v = 120  # Vehicle speed in kilometers per hour

        # Calculate realistic ranges for suspension and tire stiffness based on typical deflections
        self.mink1 = (self.m1 * 9.81) / (3.0 * 0.0254)  # Minimum suspension stiffness (spring extends about 3 inches)
        self.maxk1 = (self.m1 * 9.81) / (6.0 * 0.0254)  # Maximum suspension stiffness (softer spring)
        self.mink2 = (self.m2 * 9.81) / (1.5 * 0.0254)  # Minimum tire stiffness
        self.maxk2 = (self.m2 * 9.81) / (0.75 * 0.0254)  # Maximum tire stiffness (hard tire)

        # These hold additional analysis results
        self.accel = None  # Acceleration data (vertical)
        self.accelMax = 0.0  # Maximum vertical acceleration (in g’s)
        self.accelLim = 2.0  # Target acceleration limit (g's)
        self.SSE = 0.0  # Sum of squared error between car position and road profile


class CarView():
    def __init__(self, args):
        self.input_widgets, self.display_widgets = args
        # unpack widgets with same names as they have on the GUI
        self.le_m1, self.le_v, self.le_k1, self.le_c1, self.le_m2, self.le_k2, self.le_ang, \
         self.le_tmax, self.chk_IncludeAccel = self.input_widgets

        self.gv_Schematic, self.chk_LogX, self.chk_LogY, self.chk_LogAccel, \
        self.chk_ShowAccel, self.lbl_MaxMinInfo, self.layout_horizontal_main = self.display_widgets

        # creating a canvas to draw a figure for the car model
        self.figure = Figure(tight_layout=True, frameon=True, facecolor='none')
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.layout_horizontal_main.addWidget(self.canvas)

        # axes for the plotting using view
        self.ax = self.figure.add_subplot()
        if self.ax is not None:
            self.ax1 = self.ax.twinx()

        self.buildScene()

    def updateView(self, model=None):
        self.le_m1.setText("{:0.2f}".format(model.m1))
        self.le_k1.setText("{:0.2f}".format(model.k1))
        self.le_c1.setText("{:0.2f}".format(model.c1))
        self.le_m2.setText("{:0.2f}".format(model.m2))
        self.le_k2.setText("{:0.2f}".format(model.k2))
        self.le_ang.setText("{:0.2f}".format(model.yangdeg))
        self.le_tmax.setText("{:0.2f}".format(model.tmax))
        stTmp="k1_min = {:0.2f}, k1_max = {:0.2f}\nk2_min = {:0.2f}, k2_max = {:0.2f}\n".format(model.mink1, model.maxk1, model.mink2, model.maxk2)
        stTmp+="SSE = {:0.2f}".format(model.SSE)
        self.lbl_MaxMinInfo.setText(stTmp)
        self.doPlot(model)

    def buildScene(self):
        #create a scene object
        self.scene = qtw.QGraphicsScene()
        self.scene.setObjectName("MyScene")
        self.scene.setSceneRect(-200, -150, 400, 400)  # xLeft, yTop, Width, Heightt

        #set the scene for the graphics view object
        self.gv_Schematic.setScene(self.scene)
        #make some pens and brushes for my drawing
        self.setupPensAndBrushes()
        #Create the wheel and car body and add to scene

        self.Wheel = Wheel(0,50,50, pen=self.penWheel, wheelBrush=self.brushWheel, massBrush=self.brushMass, name = "Wheel")
        self.CarBody = MassBlock(0, -70, 100, 30, pen=self.penWheel, brush=self.brushMass, name="Car Body", mass=150)
        self.Wheel.addToScene(self.scene)
        self.scene.addItem(self.CarBody)

        #Add a spring between Car Body and Wheel
        self.spring = Spring(self.CarBody.x-20, self.CarBody.y + 15, self.Wheel.x, self.Wheel.y - 10)
        springPen = qtg.QPen(qtg.QColor("blue"))
        springPen.setWidth(2)
        self.spring.pen = springPen
        self.scene.addItem(self.spring)

        #Add a dashpot between Car Body and Wheel
        self.dashpot = qtw.QGraphicsLineItem(self.CarBody.x+30, self.CarBody.y + 15, self.Wheel.x+30, self.Wheel.y - 10)
        dashpotPen = qtg.QPen(qtg.QColor("green"))
        dashpotPen.setWidth(2)
        self.dashpot.setPen(dashpotPen)
        self.scene.addItem(self.dashpot)

        # Add the road surface with a ramp
        roadPen = qtg.QPen(qtg.QColor("black"))
        roadPen.setWidth(2)
        roadBrush = qtg.QBrush(qtg.QColor(150, 150, 150, 128))  # Gray fill for the road
        roadPath = qtg.QPainterPath()
        roadPath.moveTo(-200, self.Wheel.y + 75)  # Start at the left edge
        roadPath.lineTo(-50, self.Wheel.y + 50)  # Flat section before the ramp
        roadPath.lineTo(0, self.Wheel.y+45)  # Ramp up to the wheel's position
        roadPath.lineTo(200, self.Wheel.y+40)  # Flat section after the ramp
        roadItem = qtw.QGraphicsPathItem(roadPath)
        roadItem.setPen(roadPen)
        roadItem.setBrush(roadBrush)
        self.scene.addItem(roadItem)

        # Add the tire spring (k2) between the wheel and the road
        self.tireSpring = Spring(self.Wheel.x-10, self.Wheel.y + 20, self.Wheel.x, self.Wheel.y + 50)
        tireSpringPen = qtg.QPen(qtg.QColor("red"))
        tireSpringPen.setWidth(2)
        self.tireSpring.pen = tireSpringPen
        self.scene.addItem(self.tireSpring)
        bodyLabel = qtw.QGraphicsTextItem("Body")
        bodyLabel.setPos(self.CarBody.x - 20, self.CarBody.y - 35)
        self.scene.addItem(bodyLabel)

        #add labels for each

        suspensionLabel = qtw.QGraphicsTextItem("Suspension")
        suspensionLabel.setPos(self.CarBody.x +70, (self.CarBody.y + self.Wheel.y) / 2)
        self.scene.addItem(suspensionLabel)

        wheelLabel = qtw.QGraphicsTextItem("Wheel")
        wheelLabel.setPos(self.Wheel.x +60, self.Wheel.y - 10)
        self.scene.addItem(wheelLabel)

        roadLabel = qtw.QGraphicsTextItem("Road")
        roadLabel.setPos(self.Wheel.x +60, self.Wheel.y + 50)
        self.scene.addItem(roadLabel)

        # Add datum level line
        datumPen = qtg.QPen(qtg.QColor("black"))
        datumPen.setWidth(1)
        datumPen.setStyle(qtc.Qt.DashLine)
        datumLine = qtw.QGraphicsLineItem(-200, 125, 200, 125)
        datumLine.setPen(datumPen)
        self.scene.addItem(datumLine)

        datumLabel = qtw.QGraphicsTextItem("Datum Level")
        datumLabel.setPos(150, -250)
        self.scene.addItem(datumLabel)

    def setupPensAndBrushes(self):
        self.penWheel = qtg.QPen(qtg.QColor("orange"))
        self.penWheel.setWidth(1)
        self.brushWheel = qtg.QBrush(qtg.QColor.fromHsv(35,255,255, 64))
        self.brushMass = qtg.QBrush(qtg.QColor(200,200,200, 128))

    def doPlot(self, model=None):
        if model.results is None:
            return
        ax = self.ax
        ax1=self.ax1
        # plot result of odeint solver
        QTPlotting = True  # assumes we are plotting onto a QT GUI form
        if ax == None:
            ax = plt.subplot()
            ax1=ax.twinx()
            QTPlotting = False  # actually, we are just using CLI and showing the plot
        ax.clear()
        ax1.clear()
        t=model.timeData
        ycar = model.results[:,0]
        ywheel=model.results[:,2]
        accel=model.accelData

        if self.chk_LogX.isChecked():
            ax.set_xlim(0.001,model.tmax)
            ax.set_xscale('log')
        else:
            ax.set_xlim(0.0, model.tmax)
            ax.set_xscale('linear')

        if self.chk_LogY.isChecked():
            ax.set_ylim(0.0001,max(ycar.max(), ywheel.max()*1.05))
            ax.set_yscale('log')
        else:
            ax.set_ylim(0.0, max(ycar.max(), ywheel.max()*1.05))
            ax.set_yscale('linear')

        ax.plot(t, ycar, 'b-', label='Body Position')
        ax.plot(t, ywheel, 'r-', label='Wheel Position')
        if self.chk_ShowAccel.isChecked():
            ax1.plot(t, accel, 'g-', label='Body Accel')
            ax1.axhline(y=accel.max(), color='orange')  # horizontal line at accel.max()
            ax1.set_yscale('log' if self.chk_LogAccel.isChecked() else 'linear')

        # add axis labels
        ax.set_ylabel("Vertical Position (m)", fontsize='large' if QTPlotting else 'medium')
        ax.set_xlabel("time (s)", fontsize='large' if QTPlotting else 'medium')
        ax1.set_ylabel("Y'' (g)", fontsize = 'large' if QTPlotting else 'medium')
        ax.legend()

        ax.axvline(x=model.tramp)  # vertical line at tramp
        ax.axhline(y=model.ymag)  # horizontal line at ymag
        # modify the tick marks
        ax.tick_params(axis='both', which='both', direction='in', top=True,
                       labelsize='large' if QTPlotting else 'medium')  # format tick marks
        ax1.tick_params(axis='both', which='both', direction='in', right=True,
                       labelsize='large' if QTPlotting else 'medium')  # format tick marks
        # show the plot
        if QTPlotting == False:
            plt.show()
        else:
            self.canvas.draw()
        #update schematic parts here
        newBodyY = ycar[-1] * 100 - 150  # Scale to fit
        newWheelY = ywheel[-1] * 1 + 50


        self.CarBody.y = newBodyY
        self.Wheel.y = newWheelY

        self.CarBody.setTransform(qtg.QTransform().translate(self.CarBody.x, self.CarBody.y))
        self.Wheel.setTransform(qtg.QTransform().translate(self.Wheel.x, self.Wheel.y))

        # 🎯 UPDATE the spring and dashpot
        self.spring.updateEndpoints(self.CarBody.x, self.CarBody.y + 20, self.Wheel.x, self.Wheel.y - 20)
        self.spring.update()
        self.dashpot.setLine(self.CarBody.x, self.CarBody.y + 20, self.Wheel.x, self.Wheel.y - 20)
class CarController():
    def __init__(self, args):
        """
        This is the controller I am using for the quarter car model.
        """
        self.input_widgets, self.display_widgets = args
        #unpack widgets with same names as they have on the GUI
        self.le_m1, self.le_v, self.le_k1, self.le_c1, self.le_m2, self.le_k2, self.le_ang, \
         self.le_tmax, self.chk_IncludeAccel = self.input_widgets

        self.gv_Schematic, self.chk_LogX, self.chk_LogY, self.chk_LogAccel, \
        self.chk_ShowAccel, self.lbl_MaxMinInfo, self.layout_horizontal_main = self.display_widgets

        self.model = CarModel()
        self.view = CarView(args)

        self.chk_IncludeAccel=qtw.QCheckBox()

    def ode_system(self, X, t):
        # define the forcing function equation for the linear ramp
        # It takes self.tramp time to climb the ramp, so y position is
        # a linear function of time.
        if t < self.model.tramp:
            y = self.model.ymag * (t / self.model.tramp)
        else:
            y = self.model.ymag

        x1 = X[0]  # car position in vertical direction
        x1dot = X[1]  # car velocity  in vertical direction
        x2 = X[2]  # wheel position in vertical direction
        x2dot = X[3]  # wheel velocity in vertical direction

        # write the non-trivial equations in vertical direction
        x1ddot = (-self.model.k1*(x1-x2) - self.model.c1*(x1dot-x2dot)) / self.model.m1
        x2ddot = (self.model.k1*(x1-x2) + self.model.c1*(x1dot-x2dot) - self.model.k2*(x2-y)) / self.model.m2

        # return the derivatives of the input state vector
        return [x1dot, x1ddot, x2dot, x2ddot]

    def calculate(self, doCalc=True):
        """
        I will first set the basic properties of the car model and then calculate the result
        in another function doCalc.
        """
        #Step 1.  Read from the widgets
        self.model.m1 = float(self.le_m1.text())
        self.model.m2 = float(self.le_m2.text())
        self.model.c1 = float(self.le_c1.text())
        self.model.k1 = float(self.le_k1.text())
        self.model.k2 = float(self.le_k2.text())
        self.model.v = float(self.le_v.text())

        #recalculate min and max k values
        self.mink1=(self.model.m1*9.81)/(3.0*.0254)
        self.maxk1=(self.model.m1*9.81)/(6.0*.0254)
        self.mink2=(self.model.m2*9.81)/(1.5*.0254)
        self.maxk2=(self.model.m2*9.81)/(0.75*.0254)

        ymag=6.0/(12.0*3.3)   #This is the height of the ramp in m
        if ymag is not None:
            self.model.ymag = ymag
        self.model.yangdeg = float(self.le_ang.text())
        self.model.tmax = float(self.le_tmax.text())
        if(doCalc):
            self.doCalc()
        self.SSE((self.model.k1, self.model.c1, self.model.k2), optimizing=False)
        self.view.updateView(self.model)

    def setWidgets(self, w):
        self.view.setWidgets(w)
        self.chk_IncludeAccel=self.view.chk_IncludeAccel

    def doCalc(self, doPlot=True, doAccel=True):
        """
        This solves the differential equations for the quarter car model.
        :param doPlot:
        :param doAccel:
        :return:
        """
        v = 1000 * self.model.v / 3600  # convert speed to m/s from kph
        self.model.angrad = self.model.yangdeg * math.pi / 180.0  # convert angle to radians
        self.model.tramp = self.model.ymag / (math.sin(self.model.angrad) * v)  # calculate time to traverse ramp

        self.model.t=np.linspace(0,self.model.tmax,2000)
        ic = [0, 0, 0, 0]
        # run odeint solver
        self.model.results = odeint(self.ode_system, ic, self.model.t)
        if doAccel:
            self.calcAccel()
        self.model.timeData = self.model.t
        self.model.accelData = self.model.accel
        if doPlot:
            self.doPlot()

    def calcAccel(self):
        """
        Calculate the acceleration in the vertical direction using the forward difference formula.
        """
        N=len(self.model.t)
        self.model.accel=np.zeros(shape=N)
        vel=self.model.results[:,1]
        for i in range(N):
            if i==N-1:
                h = self.model.t[i] - self.model.t[i-1]
                self.model.accel[i]=(vel[i]-vel[i-1])/(9.81*h)  # backward difference of velocity
            else:
                h = self.model.t[i + 1] - self.model.t[i]
                self.model.accel[i] = (vel[i + 1] - vel[i]) / (9.81 * h)  # forward difference of velocity
            # else:
            #     self.model.accel[i]=(vel[i+1]-vel[i-1])/(9.81*2.0*h)  # central difference of velocity
        self.model.accelMax=self.model.accel.max()
        return True

    def OptimizeSuspension(self):
        """
        Step 1:  set parameters based on GUI inputs by calling self.set(doCalc=False)
        Step 2:  make an initial guess for k1, c1, k2
        Step 3:  optimize the suspension
        :return:
        """
        #Step 1:
        self.calculate(doCalc=False)
        #Step 2:
        x0= np.array([self.model.k1, self.model.c1, self.model.k2]) # create a numpy array with initial values for k1, c1, and k2
        #Step 3:
        answer= minimize(self.SSE, x0, method='Nelder-Mead') #use the Nelder-Mead method to minimize the SSE function (our objective function)
        self.model.k1, self.model.c1, self.model.k2 = answer.x
        self.view.updateView(self.model)

    def SSE(self, vals, optimizing=True):
        """
        Calculates the sum of square errors between the contour of the road and the car body.
        :param vals:
        :param optimizing:
        :return:
        """
        k1, c1, k2=vals  #unpack the new values for k1, c1, k2
        self.model.k1=k1
        self.model.c1=c1
        self.model.k2=k2
        self.doCalc(doPlot=False)  #solve the odesystem with the new values of k1, c1, k2
        SSE=0
        for i in range(len(self.model.results[:,0])):
            t=self.model.t[i]
            y=self.model.results[:,0][i]
            if t < self.model.tramp:
                ytarget = self.model.ymag * (t / self.model.tramp)
            else:
                ytarget = self.model.ymag
            SSE+=(y-ytarget)**2

        #some penalty functions if the constants are too small
        if optimizing:
            if k1<self.model.mink1 or k1>self.model.maxk1:
                SSE+=100
            if c1<10:
                SSE+=100
            if k2<self.model.mink2 or k2>self.model.maxk2:
                SSE+=100

            # I'm overlaying a gradient in the acceleration limit that scales with distance from a target squared.
            if self.model.accelMax>self.model.accelLim and self.chk_IncludeAccel.isChecked():
                # need to soften suspension
                SSE+=(self.model.accelMax-self.model.accelLim)**2
        self.model.SSE=SSE
        return SSE

    def doPlot(self):
        self.view.doPlot(self.model)
#endregion
#endregion

def main():
    QCM = CarController()
    QCM.doCalc()

if __name__ == '__main__':
    main()
