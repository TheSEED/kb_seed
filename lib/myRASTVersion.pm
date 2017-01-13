# This is a SAS component.
package myRASTVersion;
use base 'Class::Accessor';
__PACKAGE__->mk_accessors(qw(release));
sub new
{
    my($class) = @_;
    my $self = {
	release => "1.149",
	package_date => 1484343773,
	package_date_str => "Jan 13, 2017 15:42:53",
    };
    return bless $self, $class;
}
1;
