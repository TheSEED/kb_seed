# This is a SAS component.
package myRASTVersion;
use base 'Class::Accessor';
__PACKAGE__->mk_accessors(qw(release));
sub new
{
    my($class) = @_;
    my $self = {
	release => "1.154",
	package_date => 1506539811,
	package_date_str => "Sep 27, 2017 14:16:51",
    };
    return bless $self, $class;
}
1;
