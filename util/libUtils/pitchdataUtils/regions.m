function  [sensory,motor_premotor,STG,larynx]= regions(subject)
%Group the electrodes based on their regions from the electrode maps
%Subject   = 1 EC56
%          = 2 EC61
%          = 3 EC69
%          = 4 EC77

switch subject
    case 1
        
        motor_premotor   = sort([6:10,22:26,38:42,53:62,69:79,85:96,126:128,142:144,159:160,175:176]);
        sensory          = sort([101:107,118:125,134:141,153:158,170:174,187:192,203:208,220:224,236:240,254:256]);
        STG              = sort([1:5,17:21,33:37,49:52,65:68,81:84,97:100,113:117,129:133,145:149,161:165,177:181,193:197,209:213,225:229]);
        larynx           =  sort([motor_premotor,sensory]) ; 
    case 2        
        motor_premotor   = sort([188:192,170:176,153:157,135:141,119:124,103]);
        sensory          = sort([253:256,237:240,221:224,200:207,182:187,166:169,150,151]);
        STG              = sort([1,2,17,18,33:35,50,51,66,67,68,82:85,99:101,115:117,131:133,147:149,163:165,179:181,196,197,212:214,228:230,244:247]);
        larynx           =  sort([motor_premotor,sensory]) ; 
    case 3
        sensory   =  sort([36:45,52:64,69:79,85:93,101:103,117,118,133,134,149,150,166]);
        STG       =  sort([1:6,17:22,33:38,49:53,65:80,81:86,97:103,113:119,129:135,145:150,161:166,177:183.193:198,209:214,225:230,241:246]);
        motor_premotor         =  sort(setdiff([1:256],[sensory,STG]));
        larynx    =  sort([motor_premotor,sensory]) ; 
        
    case 4        
        sensory          =  sort([118:120,134:140,150:156,166:173,182:189,198:205,214:221,230:237,246:253]);
        STG              =  sort([1:4,17:20,33:37,49:53,65:69,81:85,97:101,113:117,129:133,145:149,161:165,177:181,193:197,209:213,225:229,241:245]);        
        motor_premotor   =  sort(setdiff([1:256],[sensory,STG]));
        larynx           =   sort([motor_premotor,sensory]) ; 
end