# BioPython-Getopt-Biological-Command-Line-Interpreter
BioPython Project using getopt library in python

![image](https://user-images.githubusercontent.com/112272836/215086629-0e56dca4-2a3d-4948-a7b4-20a23c403592.png)
![image](https://user-images.githubusercontent.com/112272836/215086685-3c20e1c6-269e-4502-9093-50e02976b2f7.png)
![image](https://user-images.githubusercontent.com/112272836/215086717-b583e80e-2ced-4742-a4aa-4d7323159985.png)
![image](https://user-images.githubusercontent.com/112272836/215086801-8fd3b15e-7730-4ac2-9c74-b21ed83713b6.png)
![image](https://user-images.githubusercontent.com/112272836/215086832-f7f6bfb1-d015-47a7-8879-223383d81fa8.png)
![image](https://user-images.githubusercontent.com/112272836/215086866-cb223302-ef46-4948-b4b6-489d6da38d01.png)
![image](https://user-images.githubusercontent.com/112272836/215086898-860cf3b1-0d39-4efb-8aef-4e899c182e39.png)
![image](https://user-images.githubusercontent.com/112272836/215086955-643446eb-df02-4c34-8715-03f975c896a7.png)
![image](https://user-images.githubusercontent.com/112272836/215086992-e541f6a0-e6c3-495c-ac29-789bbfde5612.png)
![image](https://user-images.githubusercontent.com/112272836/215087019-28b1357e-f165-4729-bc03-4c9c185a3b62.png)
![image](https://user-images.githubusercontent.com/112272836/215087060-f06b7641-a75b-46e7-b4b9-0e1908ec601a.png)
![image](https://user-images.githubusercontent.com/112272836/215087096-914e936e-c30a-46bb-9c90-08f055111491.png)

-All files used in the python code to check the functionality are uploaded including input and output results files.

-Some important remarks regarding the getopt library:

-An option followed by a colon only means that it needs an argument.

-It doesn't mean that the option is enforced. You should write your own code to enforce the existence of options/arguments.

-getopt.getopt(args, shortopts, longopts=[])

-Parses command line options and parameter list. 

-args is the argument list to be parsed, without the leading reference to the running program. Typically, this means sys.argv[1:]. 

-shortopts is the string of option letters that the script wants to recognize, 

-with options that require an argument followed by a colon (':'; i.e., the same format that Unix getopt() uses).

-exception getopt.GetoptError

-This is raised when an unrecognized option is found in the argument list or when an option requiring an argument is given none. The argument to the exception 

-is a string indicating the cause of the error. For long options, an argument given to an option which does not require one will also cause this exception to be raised. 

-The attributes msg and opt give the error message and related option; if there is no specific option to which the exception relates, opt is an empty string.

-One of the most important remarks is:

-that the getopt function reads options and forms a dictionary of each option and its argument if and only if the first 

-input was an option ,rather than that it will consider the whole input as an argument even if there are option passed and the opts array will be empty in that case.
 
-Example Usage

-Short Form

-python3 code.py -s 2000-1-2 -e 2002-1-1

-Long Form

-python3 code.py --start_date 2000-1-2 --end_date 2002-1-1

-Input Formats:

-In this code many types of inputs are allowed and also many types of errors are well handeled

-Permitted Input Formats:

-command name param1 param2 param3 ... -o option1 -x option2 -y option3 ...

-command name -o option1 -x option2 -y option3 ...

-command name param1 param2 param3 ...

-So, your input can be only options associated with their arguments parameters 

-or Only arguments

-Or a mix of options and arguments parameters

-In all these cases, many other logical errors are handeled according to the function type and functionality

-So as to give user freedom in their input format

-Also here are 20+ examples with all possible inputs from the user and using all function as well:

![1](https://user-images.githubusercontent.com/112272836/215087531-4b01a612-446d-40e2-bceb-ddd92650eb5a.PNG)
![2](https://user-images.githubusercontent.com/112272836/215087537-9f153ae9-ac2f-4cc5-8066-9c8cd46b7a96.PNG)
![3](https://user-images.githubusercontent.com/112272836/215087546-8f67424e-6ac8-48fe-9745-ed442b950c5a.PNG)
![4](https://user-images.githubusercontent.com/112272836/215087561-b13622dc-785e-4dd4-a771-17c0c43d2da8.PNG)
![5](https://user-images.githubusercontent.com/112272836/215087566-80b6da1f-822f-47cc-993b-ea20263bae10.PNG)
![6](https://user-images.githubusercontent.com/112272836/215087579-339c28ba-2c0f-40cf-b31d-3c7263ee94f1.PNG)
![7](https://user-images.githubusercontent.com/112272836/215087586-064192b0-c59e-4350-9936-c8559cba79d9.PNG)
![8](https://user-images.githubusercontent.com/112272836/215087421-af7fbadc-1ade-4dda-9931-8b545be7a89d.PNG)
![9](https://user-images.githubusercontent.com/112272836/215087430-8d8cb24d-ff1a-4a42-9e9a-54f2019fc0e4.PNG)
![10](https://user-images.githubusercontent.com/112272836/215087439-8d669541-9604-407f-baff-3b05dd0d3b65.PNG)
![11](https://user-images.githubusercontent.com/112272836/215087442-df88f59b-6d07-4632-8bb4-62422b7e5429.PNG)
![12](https://user-images.githubusercontent.com/112272836/215087444-db85e859-5062-46e5-b944-c48b08b0337b.PNG)
![13](https://user-images.githubusercontent.com/112272836/215087456-8800029e-93c6-4c45-a66f-fe8fb1e7a192.PNG)
![14](https://user-images.githubusercontent.com/112272836/215087460-43d5f1f7-44c1-4f33-a229-43af2c76d126.PNG)
![15](https://user-images.githubusercontent.com/112272836/215087464-6c1d2794-c3fd-4e71-9cc0-0313d4abba9f.PNG)
![16](https://user-images.githubusercontent.com/112272836/215087468-ddee458f-bc9b-4471-bac4-719d3d5d9dd9.PNG)
![17](https://user-images.githubusercontent.com/112272836/215087471-583ed492-2a08-43f6-a10b-082aae3ffd1a.PNG)
![18](https://user-images.githubusercontent.com/112272836/215087475-29f855e8-97c3-4cbb-80b9-2d6229a904bb.PNG)
![19](https://user-images.githubusercontent.com/112272836/215087482-e5550dc4-f9f1-4b87-b5bb-b44019a6007c.PNG)
![20](https://user-images.githubusercontent.com/112272836/215087483-9166042b-d211-49a1-b7cb-1f61af873f91.PNG)
![21](https://user-images.githubusercontent.com/112272836/215087498-64d88885-c557-4208-a85d-d15945b426eb.PNG)
![22](https://user-images.githubusercontent.com/112272836/215087510-8ecf2a93-f85e-425e-a530-77b847bf9b24.PNG)
![23](https://user-images.githubusercontent.com/112272836/215087516-21b1911b-6155-4248-b817-ec1b56d58041.PNG)
![24](https://user-images.githubusercontent.com/112272836/215087521-903dff99-1288-48ad-b821-0744cf730a41.PNG)
![25](https://user-images.githubusercontent.com/112272836/215087528-8fe416dd-67b1-44ab-8b66-0402ba2d6c3e.PNG)
