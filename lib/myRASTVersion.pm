# This is a SAS component.
package myRASTVersion;
use base 'Class::Accessor';
__PACKAGE__->mk_accessors(qw(release));
sub new
{
    my($class) = @_;
    my $self = {
	release => "1.142",
	package_date => 1480522540,
	package_date_str => "Nov 30, 2016 10:15:40",
    };
    return bless $self, $class;
}
1;
