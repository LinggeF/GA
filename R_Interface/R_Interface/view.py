from django.http import HttpResponse
from django.shortcuts import render,render_to_response
import uuid
import os
def index(request):
    return render_to_response("GA_web.html")

def upload_file(request):
    if request.method=="POST":
        #get file from request.FILES
        file=request.FILES.get("file",None)
        file_name=""
        print file
        if file == None:
            return HttpResponse(1)

        file_type = file.name.split(".")[-1]
        file_path = "upload/"
        file_name = str(uuid.uuid1()) + "." + file_type
        of = open(file_path + file_name, 'wb+')
        for chunk in file.chunks():
            of.write(chunk)
        of.close()
        return HttpResponse(file_name)

def submit_job(request):
    if request.method=="POST":
        #get parameters from request.POST
        print request.POST
        upload_file_name=request.POST.get('file_name',None)
        n_modules = request.POST.get('n_modules', None)
        ga_times = request.POST.get('ga_times', None)
        ls_times = request.POST.get('ls_times', None)
        pop_size = request.POST.get('pop_size',None)
        cros_prob = request.POST.get('cros_prob',None)
        put_prob = request.POST.get('put_prob',None)
        ri_prob = request.POST.get('ri_prob',None)
        # run r script by system call
        r_script_name="Perform_GA.R"
        base_path=os.path.abspath(".")
        upload_path=os.path.join(base_path,"upload",upload_file_name)
        output_file_name="output.txt"
        output_path=os.path.join(base_path,"output",output_file_name)

        #the format of Command Line
        cmd="Rscript"+" "+r_script_name+" "+pop_size+" "+cros_prob+" "+put_prob+" "+ga_times+" "+n_modules+" "+ri_prob+" "+ls_times+" "+upload_path+" "+output_path
        code=os.system(cmd)
        print code
        if code==0:
            return_data=0
        else:
            return_data=1
        return HttpResponse(return_data)

def download(request):
    # download
    output_path = os.path.join(os.path.abspath("."),"output","output.txt")
    with open(output_path) as f:
        c = f.read()
    return HttpResponse(c)

def hello(request):
    return HttpResponse("hello")
