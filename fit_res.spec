Name:         fit_res
Version:      1.1
Release:      alt1

Summary:      Command-line tool for fitting resonance data
Group:        System
URL:          https://github.com/slazav/fit_res
License:      GPLv3

Packager:     Vladislav Zavjalov <slazav@altlinux.org>

Source:       %name-%version.tar
BuildRequires: libgsl-devel

%description
fit_res - a command-line tool for fitting resonance data

%prep
%setup -q

%build
%make

%install
%makeinstall

%files
%_bindir/*

%changelog
* Tue Nov 14 2023 Vladislav Zavjalov <slazav@altlinux.org> 1.1-alt1
v1.1
- improve initial conditions
- shift/scale data to have values of the order of 1 in the fit
- fix error in data read, skip bad lines
- 6- or 8-parameter fit; command-line options
- --fmt_out, --do_fit, --coord options
- --pars 10 option for double-resonance fits
- update examples

* Wed Apr 07 2021 Vladislav Zavjalov <slazav@altlinux.org> 1.0-alt1
- v1.0, works for fitting 6-parameter resonances:
  (X + i*Y) = (A + i*B) + (C + i*D)/(w0^2 - w^2 - i*w*dw)
