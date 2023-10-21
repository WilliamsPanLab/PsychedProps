subjects = {'sub-MDMA001', 'sub-MDMA002', 'sub-MDMA003', 'sub-MDMA005', 'sub-MDMA007', 'sub-MDMA008', 'sub-MDMA009', 'sub-MDMA011', 'sub-MDMA012', 'sub-MDMA013', 'sub-MDMA014', 'sub-MDMA015', 'sub-MDMA016', 'sub-MDMA017'};
sessions = {'ses-00', 'ses-01', 'ses-02', 'ses-03'};

% Loop through subjects
for subjectCell = subjects
    % Loop through sessions
    for sessionCell = sessions
        % Extract the char from the cell
        subject = subjectCell{1};
        session = sessionCell{1};
        
	% Call the function with Extract_RelativeAngles(subject, session, 'rs1')
	Extract_RelativeAngles(subject, session, 'rs1');
	Extract_RelativeAngles(subject, session, 'rs2');
    end
end
