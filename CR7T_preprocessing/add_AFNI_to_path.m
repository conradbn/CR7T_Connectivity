setenv('DYLD_FALLBACK_LIBRARY_PATH',[ '/Users/conradbn/abin:' getenv('DYLD_FALLBACK_LIBRARY_PATH') ])

% is the path to AFNI already set?
[afninotfound, unixoutput] = unix('3dcopy');
% if not, then we will look for it...
if afninotfound
    % is AFNI in /abin?
    [afninotfound, unixoutput] = unix('/abin/3dcopy');
    if ~afninotfound % found it!
        setenv('PATH', [getenv('PATH') ':/abin']);
        disp(sprintf('Adding /abin to Matlab''s path (for AFNI).'));
    end
end
if afninotfound
    % is AFNI in ~/abin?
    [afninotfound, unixoutput] = unix('~/abin/3dcopy');
    if ~afninotfound % found it!
        current_directory = pwd;
        eval('cd ~');
        afni_directory = [pwd '/abin'];
        setenv('PATH', [getenv('PATH') ':' afni_directory]);
        disp(sprintf(['Adding ' afni_directory ' to Matlab''s path (for AFNI).']));
        eval(['cd ' current_directory]);
        clear current_directory afni_directory;
    end
end
if afninotfound
    % is AFNI in ~/Documents/abin?
    [afninotfound, unixoutput] = unix('~/Documents/abin/3dcopy');
    if ~afninotfound % found it!
        current_directory = pwd;
        eval('cd ~');
        afni_directory = [pwd '/Documents/abin'];
        setenv('PATH', [getenv('PATH') ':' afni_directory]);
        disp(sprintf(['Adding ' afni_directory ' to Matlab''s path (for AFNI).']));
        eval(['cd ' current_directory]);
        clear current_directory afni_directory;
    end
end
if afninotfound
    % if AFNI is not in /abin nor ~/abin, then the user will have to add
    % the correct path to AFNI before this function will work.
    disp(sprintf('*****'));
    disp(sprintf('The required AFNI functions cannot be found.'));
    disp(sprintf('AFNI needs to be added to Matlab''s path.'));
    disp(sprintf('You can do this with the following command:'));
    disp(sprintf(' '));
    disp(sprintf('   setenv(''PATH'', [getenv(''PATH'') '':/abin'']);'));
    disp(sprintf(' '));
    disp(sprintf('where (as an example) ''/abin'' is the AFNI directory.'));
    disp(sprintf('Exiting...'));
    disp(sprintf('*****'));
    return;
end
clear afninotfound unixoutput;