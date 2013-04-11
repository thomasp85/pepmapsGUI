# TODO: Add comment
# 
# Author: Thomas
###############################################################################

library(RGtk2)
test <- gClass('TestWidget', parent='GtkFrame',
		.props=list(
				gParamSpec(type='character', name='Test', nick='t', blurb='This is a test widget', default.value='Testing 1, 2')
			),
		.public=list(
				working=function(self) return('This is a public function')
			),
		.initialize=function(TestWidget){
			frame <- gtkFrame(label='test')
			frame$setBorderWidth(5)
			button <- gtkButton('test button')
			frame$add(button)
			gSignalConnect(button, 'clicked', f=function(...) print('Clicked'))
			frame$show()
			TestWidget <- frame
		}
	)


testWindow <- gtkWindow(show=FALSE)
testScroll <- gtkScrolledWindow()
testWindow$add(testScroll)
testDraw <- gtkDrawingArea()
testDraw$setEvents(794)
testScroll$addWithViewport(testDraw)
testScroll$getChild()$setSizeRequest(width=500, height=500)
testWindow$show()
asCairoDevice(testDraw)
plot(1:10, 1:10)
gSignalConnect(testDraw, 'expose-event', f=function(widget, event){
			asCairoDevice(widget)
			widget$setData(key='Device', data=dev.cur())
			return(FALSE)
		})
id <- gSignalConnect(testDraw, 'motion-notify-event', f=function(widget, event){
			if(event$state == 1024){
				pos <- widget$getPointer()
				diff <- sqrt((event$x-pos$x)^2+(event$y-pos$y)^2)
				oldSize <- widget$getAllocation()$allocation
				absSize <- 10*diff/sqrt(oldSize$width^2+oldSize$height^2)
				if(event$y-pos$y > 0){
					newSize <- lapply(oldSize, function(x) x+x*absSize)
				} else {
					newSize <- lapply(oldSize, function(x) x-x*absSize)
				}
				widget$setSizeRequest(height=newSize$height, width=newSize$width)
			} else if(event$state == 256){
				pos <- widget$getPointer()
				viewport <- widget$getParent()
				#viewportPlotRatio <- c(x=viewport$getAllocation()$allocation$width/widget$getSizeRequest()$width, y=viewport$getAllocation()$allocation$height/widget$getSizeRequest()$height)
				hChange <- viewport$getHadjustment()$getValue()-(pos$x-event$x)
				vChange <- viewport$getVadjustment()$getValue()-(pos$y-event$y)
				if(hChange >= viewport$getHadjustment()$getLower() & hChange <= viewport$getHadjustment()$getUpper()-viewport$getAllocation()$allocation$width){
					viewport$getHadjustment()$setValue(hChange)
					viewport$getHadjustment()$valueChanged()
				} else {}
				if(vChange >= viewport$getVadjustment()$getLower() & vChange <= viewport$getVadjustment()$getUpper()-viewport$getAllocation()$allocation$height){
					viewport$getVadjustment()$setValue(vChange)
					viewport$getVadjustment()$valueChanged()
				} else {}
			}
			return(FALSE)
		})

testWindow$show()