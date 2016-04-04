function capture(cmd, varargin)

    persistent prefix count;
    
    if nargin < 1, cmd = ''; end;
    
    switch cmd
        
        case 'setup'
            prefix = varargin{1};
            count  = 0;
            
        case ''
            count = count + 1;
            print([prefix '_' num2str(count) '.png'], '-dpng');
            
        otherwise
            error('?');
            
    end
    
end
