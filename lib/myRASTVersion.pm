# This is a SAS component.
package myRASTVersion;
use base 'Class::Accessor';
__PACKAGE__->mk_accessors(qw(release));
sub new
{
    my($class) = @_;
    my $self = {
	release => "1.159",
	package_date => 1506546050,
	package_date_str => "Sep 27, 2017 16:00:50",
    };
    return bless $self, $class;
}
1;
