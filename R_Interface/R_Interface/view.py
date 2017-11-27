from django.http import HttpResponse
from django.shortcuts import render,render_to_response
import uuid
import os
import json


def index(request):
    return render_to_response("GA_web.html")
def handle_upload_file(file,forlder_path):
    of = open(os.path.join(forlder_path,file.name), 'wb+')
    for chunk in file.chunks():
        of.write(chunk)
    of.close()

def upload_file(request):
    if request.method=="POST":
        #get file from request.FILES
        files=[]
        f_list=request.FILES
        for i in range(0,len(f_list)):
            files.append(f_list.get(str(i)))
        print f_list.get(str(0))
        print files
        if files == None:
            return HttpResponse(1)
        base_path=os.path.abspath(".")
        u_path = str(uuid.uuid1())
        folder_path = os.path.join(base_path,"upload/", u_path)
        if not os.path.exists(folder_path):
            os.mkdir(folder_path)
        for f in files:
            handle_upload_file(f,folder_path)
        return HttpResponse(str(folder_path))

def submit_job(request):
    if request.method=="POST":
        #get parameters from request.POST
        print request.POST
        upload_folder_path=request.POST.get('folder_path',None)
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
        input1=os.path.join(upload_folder_path,"Score_matrix_even_inter.txt")
        input2=os.path.join(upload_folder_path,"mR_mR_corr.txt")
        input3=os.path.join(upload_folder_path,"gene_category.txt")
        output_file_name="output.txt"
        output_path=os.path.join(base_path,"output","output.txt")

        #the format of Command Line
        cmd="Rscript"+" "+r_script_name+" "+pop_size+" "+cros_prob+" "+put_prob+" "+ga_times+" "+n_modules+" "+ri_prob+" "+ls_times+" "+input1+" "+input2+" "+input3+" "+output_path
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
