% load in mouse DMN
networks = load('/oak/stanford/groups/leanew1/users/apines/p50_mice/Mouse_DMN_components.mat');
Dnet = networks.components_matrix';
Mask = networks.mask;

% re-boolean and transpose mask for consistency
Mask = Mask == 1;
Mask = Mask';

% treshold DMN mask (threshold > 0.6)
DMN_bool = Dnet > 0.6;
DMN_bool(Mask == 0) = 0;
DMN_bool = logical(DMN_bool);

% Apply masking rules
Dnet(~DMN_bool) = 0;        % Gray for non-DMN brain areas
Dnet(Mask == 0) = -0.1;     % White for non-brain areas


% need to rotate it around centerpoint... back edge of DMN should be on cortex-non-cortex border and project inwards
