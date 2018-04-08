classdef logging < handle
    %
    % NOSET 5
    % DEBUG 0
    % INFO 1
    % WARNING 2
    % ERROR 3
    % CRITICAL 4
    
    properties (Access = private)
        level
    end
    
    methods
        function this = logging()
            this.level = 2;
        end
        
        function setlevel(this, lvl)
            this.level = lvl;
        end
        
        function info(this, s)
            if (this.level <= 1)
                out = strcat({datestr(datetime('now'))}, {' '}, {s});
                cprintf('*[0 ,0.8, 0]', out{1});
                cprintf('\n');
            end
        end
        
        function warning(this, s) 
            if (this.level <= 2)
                out = strcat({datestr(datetime('now'))}, {' '}, {s});
                cprintf('*[1, 0.6, 0.2]', out{1});
                cprintf('\n');
            end
        end
        
        function critical(this, s)
            if (this.level <= 4)
                out = strcat({datestr(datetime('now'))}, {' '}, {s});
                cprintf('*[0, 1., 1.]',  out{1});
                cprintf('\n');
            end
        end
        
        function debug(this, s)
            if (this.level <= 0)
                out = strcat({datestr(datetime('now'))}, {' '}, {s});
                cprintf('*[1, 1, 0.2]',  out{1});
                cprintf('\n');
            end
        end
        
        function error(this, s)
            if (this.level <= 3)
                out = strcat({datestr(datetime('now'))}, {' '}, {s});
                cprintf('*[1, 0 ,0.]', out{1});
                cprintf('\n');
            end
        end
    end
    
    methods (Static)      

    end
    
end

