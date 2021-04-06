Name:         fit_res
Version:      1.0
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
* Wed Apr 07 2021 Vladislav Zavjalov <slazav@altlinux.org> 1.0-alt1
- v1.0, works for fitting 6-parameter resonances:
  (X + i*Y) = (A + i*B) + (C + i*D)/(w0^2 - w^2 - i*w*dw)
