Summary: Gosia coulex cross-section code
%global version 20110524.4
%global release 1
Name: gosia
Version: %{version}
Release: %{release}%{dist}
License: public domain
Vendor: D. cline, T. Czosnyka, A.B. Hayes, P. Napiorkowski, N. Warr, C.Y. Wu
Group: Applications/Analysis
Source: /usr/src/redhat/SOURCES/gosia.tar.gz
Prefix:/usr
BuildRoot: /tmp/package_%{name}-%{version}.%{release}
BuildRequires: gcc-gfortran

%description
Gosia is a code for calculating cross-sections for Coulomb excitation
experiments.

%prep
%setup

%build
make %{?_smp_mflags}

%install
%make_install

%files
%defattr(-,root,root,-)
/usr/bin/*
/usr/share/man/man1/*
