import webapp2 

def applyChangesVideo():	#applies changes to video/iamge, returns changed video
	return None

def score(fps, px_size):	#calculates movie score and returns integer score
	return 100; 

def redflags(): # returns string of redflags 
	return "warning: motion blur"

def displayVideo(some video file): #takes in a video file and displays it on the page. Use for both input 
                    #and output video 
    return None

class MainPage(webapp2.RequestHandler):
	def get(self):
		fps = self.request.get("fps")
		px_size = self.request.get("px_size")
		video_score = score(fps, px_size)
		red_flags = redflags()
        input_video = displayVideo()
        output_video = displayVideo(applyChangesVideo());
		self.response.headers["Content-Type"] = "text/html"
		self.response.write(""" website.html stuff """.format(input_video, fps, px_size, score, output_video, red_flags))

routes = [('/', MainPage)]

my_app = webapp2.WSGIApplication(routes, debug = True)


