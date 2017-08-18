# This is a SAS component.
package myRASTVersion;
use base 'Class::Accessor';
__PACKAGE__->mk_accessors(qw(release));
sub new
{
    my($class) = @_;
    my $self = {
	release => "1.152",
	package_date => 1503079606,
	package_date_str => "Aug 18, 2017 13:06:46",
    };
    return bless $self, $class;
}
1;
