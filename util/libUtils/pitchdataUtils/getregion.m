function eregion=getregion(subject,electrode_number)
% finds which region electrode_number of subject  belongs to
% eregion = 1 motor_premotor
%         = 2 sensory
%         = 3 STG
%         = 4 none of the above
[sensory,motor_premotor,STG,larynx]= regions(subject);
if find(motor_premotor==electrode_number)
    eregion=1;
    fprintf('\n Electrode %d belongs to %s\n',electrode_number,'motor_premotor');
elseif find(sensory==electrode_number)
    eregion=2;
    fprintf('\n Electrode %d belongs to %s\n',electrode_number,'sensory');
elseif find(STG==electrode_number)
    eregion=3;
    fprintf('\n Electrode %d belongs to %s\n',electrode_number,'STG');
else
    eregion=4;
    fprintf('\n Electrode %d doesn''t belong to any of the regions of interest')
end

end