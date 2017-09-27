# This is a SAS component.
package myRASTVersion;
use base 'Class::Accessor';
__PACKAGE__->mk_accessors(qw(release));
sub new
{
    my($class) = @_;
    my $self = {
	release => "1.158",
	package_date => 1506542994,
	package_date_str => "Sep 27, 2017 15:09:54",
    };
    return bless $self, $class;
}
1;
