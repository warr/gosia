Summary: Gosia coulex cross-section code
%define version 20090915
%define release 1
Name: gosia
Version: %{version}
Release: %{release}%{dist}
License: public domain
Vendor: Nigel Warr
Group: Applications/Analysis
Source: /usr/src/redhat/SOURCES/gosia.tar.gz
Prefix:/usr
BuildRoot: /tmp/package_%{name}-%{version}.%{release}
Requires: gosia_doc

%description
Gosia is a code for calculating cross-sections for Coulomb excitation
experiments.

%prep
%setup

%build
make

%install
make ROOT="$RPM_BUILD_ROOT" install

%files
%defattr(-,root,root)
/usr/bin/gosia
/usr/share/man/man1/gosia.1.gz

%package -n gosia_doc
Summary: Documentation for gosia, gosia2, pawel
Group: Application/Physics
%description -n gosia_doc
Documentation for gosia, gosia2, pawel etc. This package contains the
scanned-in original manual and a version typed up from that in 2007 by Adam
Hayes.

%files -n gosia_doc
%defattr(-,root,root)
%doc doc/2008manual.pdf
